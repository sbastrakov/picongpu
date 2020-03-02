/* Copyright 2013-2020 Axel Huebl, Heiko Burau, Rene Widera, Richard Pausch,
 *                     Klaus Steiniger, Felix Schmitt, Benjamin Worpitz,
 *                     Juncheng E, Sergei Bastrakov
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "picongpu/simulation_defines.hpp"

#include "picongpu/particles/traits/SpeciesEligibleForSolver.hpp"
#include "picongpu/plugins/ISimulationPlugin.hpp"
#include "picongpu/plugins/common/stringHelpers.hpp"
#include "picongpu/plugins/saxs/Saxs.kernel"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/dimensions/DataSpaceOperations.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>
#include <pmacc/memory/MakeUnique.hpp>
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/mpi/reduceMethods/Reduce.hpp>
#include <pmacc/nvidia/functors/Add.hpp>
#include <pmacc/traits/GetNumWorkers.hpp>
#include <pmacc/traits/HasIdentifier.hpp>

#include <boost/filesystem.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace picongpu
{
namespace plugins
{
namespace saxs
{
namespace detail
{

    /** Writer of X-ray diffraction results to files
     *
     * Intended to be used from a single rank in a simulation
     */
    class Writer
    {
    public:

        Writer( std::string const & prefix ):
            prefix( prefix ),
            directoryName( "xrayDiffraction" )
        {
            auto & fs = Environment< simDim >::get().Filesystem();
            fs.createDirectory( directoryName );
            fs.setDirectoryPermissions( directoryName );
        }

        void write(
            ReciprocalSpace const & reciprocalSpace,
            std::vector< float_X > const & intensity,
            float_X const combinedWeighting,
            int64_t const numMacroparticles,
            uint32_t const currentStep
        ) const
        {
            auto const baseFileName = directoryName + prefix +
                "_" + std::to_string( currentStep );
            writeIntensity(
                reciprocalSpace,
                intensity,
                baseFileName + ".dat"
            );
            writeLog(
                combinedWeighting,
                numMacroparticles,
                baseFileName + ".log"
            );
        }

    private:

        //! Directory name inside the general output directory
        std::string const directoryName;

        std::string const prefix;

        /**
        * Write scattering intensity for each q
        * @param intensity
        * @param name The name of output file
        **/
        void writeIntensity(
            ReciprocalSpace const & reciprocalSpace,
            std::vector< float_X > const & intensity,
            std::string fileName
        ) const
        {
            auto ofile = std::ofstream{ fileName.c_str() };
            if (!ofile)
            {
                std::cerr << "Can't open file [" << fileName
                    << "] for output, disable plugin output.\n";
            }
            else
            {
                ofile << intensity.size() << "\n";
                ofile << "# qx qy qz intensity \n";
                for( int i = 0; i < intensity.size(); i++ )
                    ofile << reciprocalSpace.getValue( i ) << " " << intensity[ i ] << "\n";
                ofile.close();
            }
        }

        /**
        * Write a log file for number of real particles and number of
        * macro particles.
        * @param name The name of output file
        **/
        void writeLog( 
            float_X const combinedWeighting,
            int64_t const numMacroparticles,
            std::string & const fileName
        ) const
        {
            auto ofile = std::ofstream{ fileName.c_str() };
            if( !ofile )
            {
                std::cerr << "Can't open file [" << fileName
                    << "] for X-ray diffraction output, disable plugin output.\n";
            }
            else
            {
                ofile << "Number of particles (combined weighting):"
                    << " " << combinedWeighting << "\n";
                ofile << "Number of macroparticles:"
                    << " " << numMacroparticles << "\n";
                ofile.close();
            }
        }

    };

    template< typename T_Species >
    class Implementation
    {
    public:

        Implementation(
            ReciprocalSpace & const reciprocalSpace,
            std::string & const prefix
        ):
            reciprocalSpace( reciprocalSpace ),
            prefix( prefix )
        {
            isMasterRank = reduce.hasResult(mpi::reduceMethods::Reduce());
            auto totalNumVectors = reciprocalSpace.size.productOfComponents( );
            auto size = DataSpace< DIM1 >( totalNumVectors );
            sumfcoskr = pmacc::memory::makeUnique< FloatBuffer >( size );
            sumfsinkr = pmacc::memory::makeUnique< FloatBuffer >( size );
            // allocate one float on GPU and host
            combinedWeighting = pmacc::memory::makeUnique< FloatBuffer >( DataSpace< DIM1 >( 1 ) );
            // allocate one int on GPU and host
            numMacroparticles = pmacc::memory::makeUnique< IntBuffer >( DataSpace< DIM1 >( 1 ) );
            globalSumfcoskr.resize( totalNumVectors );
            globalSumfsinkr.resize( totalNumVectors );
            intensity.resize( totalNumVectors );

            // only rank 0 create a file
            if( isMasterRank )
            {
                pmacc::Filesystem< simDim > & fs = Environment<simDim>::get().Filesystem();
                fs.createDirectory( "saxsOutput" );
                fs.setDirectoryPermissions( "saxsOutput" );
                writer = pmacc::memory::makeUnique< Writer >( prefix );
            }

        }

        void run(
            uint32_t const currentStep,
            MappingDesc & const cellDescription
        )
        {
            computeDiffraction( cellDescription );
            reduceResults();
            if( isMasterRank )
            {
                computeIntensity();
                writer->write(
                    reciprocalSpace,
                    intensity,
                    globalCombinedWeighting,
                    globalCombinedNumMacroparticles,
                    currentStep
                );
            }
        }

    private:

        using FloatBuffer = GridBuffer<
            float_X,
            DIM1
        >;
        using IntBuffer = GridBuffer<
            int64_t,
            DIM1
        >;

        //! Results for local domain

        //! The real part of structure factor
        std::unique_ptr< FloatBuffer > sumfcoskr;
        //! The imaginary part of structure factor
        std::unique_ptr< FloatBuffer > sumfsinkr;
        //! Number of real particles
        std::unique_ptr< FloatBuffer > combinedWeighting;
        //! Number of macro particles
        std::unique_ptr< IntBuffer > numMacroparticles;

        ReciprocalSpace reciprocalSpace;

        //! Reduced results for the global domain
        std::vector< float_X > globalSumfcoskr;
        std::vector< float_X > globalSumfsinkr;
        std::vector< float_X > intensity;
        float_X globalCombinedWeighting;
        int64_t globalCombinedNumMacroparticles;

        std::unique_ptr< Writer > writer;

        bool isMasterRank;

        mpi::MPIReduce reduce;
        std::string prefix;

        /** Compute diffration for macroparticles of the local domain
         *
         * @param cellDescription mapping description
         */
        void computeDiffraction( MappingDesc & const cellDescription )
        {
            // calculate the absolute position of the particles
            auto const & subGrid = Environment< simDim >::get().SubGrid();
            auto const localDomainOffset = subGrid.getLocalDomain().offset;

            constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                pmacc::math::CT::volume< MappingDesc::SuperCellSize >::type::value
            >::value;

            // initialize variables with zero
            sumfcoskr->getDeviceBuffer( ).setValue( 0.0 );
            sumfsinkr->getDeviceBuffer( ).setValue( 0.0 );
            combinedWeighting->getDeviceBuffer( ).setValue( 0.0 );
            numMacroparticles->getDeviceBuffer( ).setValue( 0.0 );          

            // PIC-like kernel call of the SAXS kernel
            DataConnector &dc = Environment<>::get().DataConnector();
            auto particles =
                dc.get<T_Species>(T_Species::FrameType::getName(), true);
            auto const totalNumVectors = reciprocalSpace.size.productOfComponents( );
            auto const numBlocks = ( totalNumVectors + numWorkers - 1 ) / numWorkers;
            PMACC_KERNEL(
                KernelSaxs< numWorkers >{ }
            )(
                numBlocks,
                numWorkers
                )(
                    // Pointer to particles memory on the device
                    particles->getDeviceParticlesBox( ),
                    // Pointer to memory of sumfcoskr & sumfsinkr on the device
                    sumfcoskr->getDeviceBuffer( ).getDataBox( ),
                    sumfsinkr->getDeviceBuffer( ).getDataBox( ),
                    combinedWeighting->getDeviceBuffer( ).getDataBox( ),
                    numMacroparticles->getDeviceBuffer( ).getDataBox( ),
                    localDomainOffset,
                    cellDescription,
                    reciprocalSpace
                    );

            dc.releaseData( T_Species::FrameType::getName( ) );

            sumfcoskr->deviceToHost();
            sumfsinkr->deviceToHost();
            combinedWeighting->deviceToHost();
            numMacroparticles->deviceToHost();
            __getTransactionEvent().waitForFinished();
        }

        /** Collect intensity data from each CPU and store result on master
        *  copyIntensityDeviceToHost should be called before */
        void reduceResults()
        {
            auto const totalNumVectors = reciprocalSpace.size.productOfComponents( );
            globalCombinedWeighting = 0._X;
            globalCombinedNumMacroparticles = 0;
            reduce(
                nvidia::functors::Add( ),
                globalSumfcoskr.data( ),
                sumfcoskr->getHostBuffer( ).getBasePointer( ),
                totalNumVectors,
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                globalSumfsinkr.data( ),
                sumfsinkr->getHostBuffer( ).getBasePointer(),
                totalNumVectors,
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                &globalCombinedWeighting,
                combinedWeighting->getHostBuffer( ).getBasePointer( ),
                1,
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                &globalCombinedNumMacroparticles,
                numMacroparticles->getHostBuffer( ).getBasePointer( ),
                1,
                mpi::reduceMethods::Reduce( )
            );
        }

        void computeIntensity()
        {
            auto const size = intensity.size();
            for( int i = 0; i < size; i++ )
                intensity[i] =
                (globalSumfcoskr[i] * globalSumfcoskr[i] +
                    globalSumfsinkr[i] * globalSumfsinkr[i]) /
                globalCombinedWeighting;
        }


 

    };

} // namespace detail

    using namespace pmacc;

    namespace po = boost::program_options;

    /** SAXS plugin
     * This SAXS plugin simulates the SAXS scattering pattern on-the-fly
     * from the particle positions obtained from PIConGPU.
     **/
    /// XrayDiffraction
    template< typename ParticlesType >
    class Saxs : public ISimulationPlugin
    {
    private:

        using Implementation = detail::Implementation< ParticlesType >;
        std::unique_ptr< Implementation > pImpl;

        MappingDesc cellDescription;
        std::string notifyPeriod; 
        std::string prefix;

        /** Range of scattering vector
         * The scattering vector here is defined as
         * 4*pi*sin(theta)/lambda, where 2theta is the angle
         * between scattered and incident beam
         **/
        float3_X q_min, q_max;

        //! Number of scattering vectors
        DataSpace< 3 > numVectors;

    public:

        Saxs() :           
            prefix( ParticlesType::FrameType::getName() + std::string("_saxs") ),
            cellDescription( DataSpace< simDim >::create( 0 ) )
        {
            Environment<>::get().PluginConnector().registerPlugin( this );
        }

        /**
         * This function represents what is actually calculated if the plugin
         * is called. Here, one only sets the particles pointer to the data of
         * the latest time step and calls the 'calculateSAXS'
         * function if for the actual time step radiation is to be calculated.
         * @param currentStep
         */
        void notify(uint32_t currentStep)
        {
            pImpl->run(
                currentStep,
                cellDescription
            );    
        }

        void pluginRegisterHelp( po::options_description &desc ) override
        {
            desc.add_options()((prefix + ".period").c_str(),
                po::value<std::string>(&notifyPeriod),
                "enable plugin [for each n-th step]")(
                (prefix + ".qx_max").c_str(),
                po::value<float_X>(&q_max[0])->default_value(5),
                "reciprocal space range qx_max (A^-1)")(
                (prefix + ".qy_max").c_str(),
                po::value<float_X>(&q_max[1])->default_value(5),
                "reciprocal space range qy_max (A^-1)")(
                (prefix + ".qz_max").c_str(),
                po::value<float_X>(&q_max[2])->default_value(5),
                "reciprocal space range qz_max (A^-1)")(
                (prefix + ".qx_min").c_str(),
                po::value<float_X>(&q_min[0])->default_value(-5),
                "reciprocal space range qx_min (A^-1)")(
                (prefix + ".qy_min").c_str(),
                po::value<float_X>(&q_min[1])->default_value(-5),
                "reciprocal space range qy_min (A^-1)")(
                (prefix + ".qz_min").c_str(),
                po::value<float_X>(&q_min[2])->default_value(-5),
                "reciprocal space range qz_min (A^-1)")(
                (prefix + ".n_qx").c_str(),
                po::value<int>(&numVectors[0])->default_value(100),
                "Number of qx")((prefix + ".n_qy").c_str(),
                po::value<int>(&numVectors[1])->default_value(100),
                "Number of qy")((prefix + ".n_qz").c_str(),
                po::value<int>(&numVectors[2])->default_value(1), "Number of qz");
        }

        std::string pluginGetName() const override
        {
            return std::string{ "SAXS: calculate the SAXS scattering intensity of a species" };
        }

        void setMappingDescription( MappingDesc * newCellDescription ) override
        {
            cellDescription = *newCellDescription;
        }

        void restart(
            uint32_t restartStep,
            std::string const restartDirectory
        ) override
        {
            // No state to be read
        }

        void checkpoint(
            uint32_t currentStep,
            std::string const restartDirectory
        ) override
        {
            // No state to be saved
        }

    private:

        /**
         * The plugin is loaded on every MPI rank, and therefore this function is
         * executed on every MPI rank.
         * One host with MPI rank 0 is defined to be the master.
         * It creates a folder where all the
         * results are saved in a plain text format.
         **/
        void pluginLoad()
        {
            if (!notifyPeriod.empty())
            {
                auto qStep = ( q_max - q_min ) / precisionCast< float_X >( numVectors );
                auto reciprocalSpace = detail::ReciprocalSpace{
                    q_min,
                    qStep,
                    numVectors
                };
                pImpl = pmacc::memory::makeUnique< Implementation >(
                    reciprocalSpace,
                    prefix
                );
                Environment<>::get().PluginConnector().setNotificationPeriod(
                    this,
                    notifyPeriod
                );
            }
        }

        void pluginUnload() override
        {
            if( !notifyPeriod.empty() )
                CUDA_CHECK( cudaGetLastError() );
        }

    };

} // namespace saxs
} // namespace plugins

namespace particles
{
namespace traits
{
    template<
        typename T_Species,
        typename T_UnspecifiedSpecies
    >
    struct SpeciesEligibleForSolver<
        T_Species,
        plugins::saxs::Saxs< T_UnspecifiedSpecies >
    >
    {
        using FrameType = typename T_Species::FrameType;

        using RequiredIdentifiers = MakeSeq_t<
            localCellIdx,
            position<>,
            weighting
        >;

        using type = typename pmacc::traits::HasIdentifiers<
            FrameType,
            RequiredIdentifiers
        >::type;
    };

} // namespace traits
} // namespace particles
} // namespace picongpu
