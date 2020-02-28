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
        using SuperCellSize = MappingDesc::SuperCellSize;

        using FloatBuffer = GridBuffer< float_X, DIM1 >;
        using IntBuffer = GridBuffer< int64_t, DIM1 >;

        //! The real part of structure factor
        std::unique_ptr< FloatBuffer > sumfcoskr;
        //! The imaginary part of structure factor
        std::unique_ptr< FloatBuffer > sumfsinkr;
        //! Number of real particles
        std::unique_ptr< FloatBuffer > np;
        //! Number of macro particles
        std::unique_ptr< IntBuffer > nmp;

        MappingDesc cellDescription;
        std::string notifyPeriod; 
        std::string prefix;

        /** Range of scattering vector
         * The scattering vector here is defined as
         * 4*pi*sin(theta)/lambda, where 2theta is the angle
         * between scattered and incident beam
         **/
        float3_X q_min, q_max;
        detail::ReciprocalSpace reciprocalSpace;

        //! Number of scattering vectors
        DataSpace< 3 > numVectors;

        std::vector< float_X > sumfcoskr_master;
        std::vector< float_X > sumfsinkr_master;
        std::vector< float_X > intensity_master;
        float_X np_master;
        int64_t nmp_master;

        bool isMaster;

        uint32_t currentStep;

        mpi::MPIReduce reduce;

    public:

        Saxs() :           
            prefix( ParticlesType::FrameType::getName() + std::string("_saxs") ),
            isMaster(false),
            currentStep(0),
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
            auto qStep = ( q_max - q_min ) / precisionCast< float_X >( numVectors );
            reciprocalSpace = detail::ReciprocalSpace{
                q_min,
                qStep,
                numVectors
            };
            calculateSAXS(currentStep);
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

        void setMappingDescription( MappingDesc * cellDescription )
        {
            this->cellDescription = *cellDescription;
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
                isMaster = reduce.hasResult(mpi::reduceMethods::Reduce());
                auto totalNumVectors = numVectors.productOfComponents( );
                auto size = DataSpace< DIM1 >( totalNumVectors );
                sumfcoskr = pmacc::memory::makeUnique< FloatBuffer >( size );
                sumfsinkr = pmacc::memory::makeUnique< FloatBuffer >( size );
                // allocate one float on GPU and host
                np = pmacc::memory::makeUnique< FloatBuffer >( DataSpace< DIM1 >( 1 ) );
                // allocate one int on GPU and host
                nmp = pmacc::memory::makeUnique< IntBuffer >( DataSpace< DIM1 >( 1 ) );
                sumfcoskr_master.resize( totalNumVectors );
                sumfsinkr_master.resize( totalNumVectors );
                intensity_master.resize( totalNumVectors );
                Environment<>::get().PluginConnector().setNotificationPeriod(
                    this, notifyPeriod);

                // only rank 0 create a file
                if (isMaster)
                {
                    pmacc::Filesystem<simDim> & fs =
                        Environment<simDim>::get().Filesystem();
                    fs.createDirectory( "saxsOutput" );
                    fs.setDirectoryPermissions( "saxsOutput" );
                }
            }
        }

        /**
         * Write scattering intensity for each q
         * @param intensity
         * @param name The name of output file
         **/
        void writeIntensity(const std::vector< float_X > & intensity, std::string name)
        {
            std::ofstream ofile;
            ofile.open(name.c_str(), std::ofstream::out | std::ostream::trunc);
            if (!ofile)
            {
                std::cerr << "Can't open file [" << name
                          << "] for output, disable plugin output.\n";
                isMaster = false; // no Master anymore -> no process is able to write
            }
            else
            {
                auto const totalNumVectors = numVectors.productOfComponents( );
                ofile << totalNumVectors << "\n";
                ofile << "# qx qy qz intensity \n";
                for( int i = 0; i < totalNumVectors; i++ )
                    ofile << reciprocalSpace.getValue( i ) << " " << intensity[ i ] << "\n";
                ofile.close();
            }
        }

        /**
         * Write a log file for number of real particles and number of
         * macro particles.
         * @param np_master The number of real particles
         * @param nmp_master The number of macro particles
         * @param name The name of output file
         **/
        void writeLog(
            float_X np_master,
            int64_t nmp_master,
            std::string name
        )
        {
            auto ofile = std::ofstream{ name.c_str() };
            if( !ofile )
            {
                std::cerr << "Can't open file [" << name
                          << "] for output, disable plugin output.\n";
                isMaster = false; // no Master anymore -> no process is able to write
            }
            else
            {
                ofile << "Number of particles:"
                      << " " << np_master << "\n";
                ofile << "Number of macro particles:"
                      << " " << nmp_master << "\n";
                ofile.close();
            }
        }

        void pluginUnload()
        {
            if (!notifyPeriod.empty())
            {
                CUDA_CHECK(cudaGetLastError());
            }
        }

        //! Method to copy data from GPU to CPU
        void copyIntensityDeviceToHost()
        {
            sumfcoskr->deviceToHost();
            sumfsinkr->deviceToHost();
            np->deviceToHost();
            nmp->deviceToHost();
            __getTransactionEvent().waitForFinished();
        }

        /** Collect intensity data from each CPU and store result on master
         *  copyIntensityDeviceToHost should be called before */
        void collectIntensityOnMaster()
        {
            auto const totalNumVectors = numVectors.productOfComponents( );
            reduce(
                nvidia::functors::Add( ),
                sumfcoskr_master.data( ),
                sumfcoskr->getHostBuffer( ).getBasePointer( ),
                totalNumVectors,
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                sumfsinkr_master.data( ),
                sumfsinkr->getHostBuffer( ).getBasePointer(),
                totalNumVectors,
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                &np_master,
                np->getHostBuffer( ).getBasePointer( ),
                1,
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                &nmp_master,
                nmp->getHostBuffer( ).getBasePointer( ),
                1,
                mpi::reduceMethods::Reduce( )
            );

            // Calculate intensity on master
            if (isMaster)
            {
                for( int i = 0; i < totalNumVectors; i++ )
                    intensity_master[i] =
                        (sumfcoskr_master[i] * sumfcoskr_master[i] +
                            sumfsinkr_master[i] * sumfsinkr_master[i]) /
                        np_master;

                std::stringstream o_step;
                o_step << currentStep;
                writeIntensity(intensity_master,
                    "saxsOutput/" + prefix + "_" + o_step.str() + ".dat");
                writeLog(np_master, nmp_master,
                    "saxsOutput/" + prefix + "_" + o_step.str() + ".log");
            }
        }

        // perform all operations to get data from GPU to master
        void collectDataGPUToMaster()
        {
            // collect data GPU -> CPU -> Master
            copyIntensityDeviceToHost();
            collectIntensityOnMaster();
        }

        /**
         * This functions calls the SAXS kernel. It specifies how the
         * calculation is parallelized.
         **/
        void calculateSAXS(uint32_t currentStep)
        {
            this->currentStep = currentStep;

            DataConnector &dc = Environment<>::get().DataConnector();
            auto particles =
                dc.get<ParticlesType>(ParticlesType::FrameType::getName(), true);

            // calculate the absolute position of the particles
            const SubGrid<simDim> &subGrid = Environment<simDim>::get().SubGrid();
            DataSpace<simDim> globalOffset(subGrid.getLocalDomain().offset);

            constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                pmacc::math::CT::volume< SuperCellSize >::type::value
            >::value;

            // initialize variables with zero
            sumfcoskr->getDeviceBuffer( ).setValue( 0.0 );
            sumfsinkr->getDeviceBuffer( ).setValue( 0.0 );
            np->getDeviceBuffer( ).setValue( 0.0 );
            nmp->getDeviceBuffer( ).setValue( 0.0 );          

            // PIC-like kernel call of the SAXS kernel
            auto const totalNumVectors = numVectors.productOfComponents( );
            auto const numBlocks = ( totalNumVectors + numWorkers - 1 ) / numWorkers;
            PMACC_KERNEL(
                detail::KernelSaxs< numWorkers >{ }
            )(
                numBlocks,
                numWorkers
            )(
                // Pointer to particles memory on the device
                particles->getDeviceParticlesBox( ),
                // Pointer to memory of sumfcoskr & sumfsinkr on the device
                sumfcoskr->getDeviceBuffer( ).getDataBox( ),
                sumfsinkr->getDeviceBuffer( ).getDataBox( ),
                np->getDeviceBuffer( ).getDataBox( ),
                nmp->getDeviceBuffer( ).getDataBox( ),
                globalOffset,
                cellDescription,
                reciprocalSpace
            );

            dc.releaseData( ParticlesType::FrameType::getName( ) );

            collectDataGPUToMaster();

            // reset amplitudes on GPU back to zero
            sumfcoskr->getDeviceBuffer( ).reset( false );
            sumfsinkr->getDeviceBuffer( ).reset( false );
            np->getDeviceBuffer( ).reset( false );
            nmp->getDeviceBuffer( ).reset( false );
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
