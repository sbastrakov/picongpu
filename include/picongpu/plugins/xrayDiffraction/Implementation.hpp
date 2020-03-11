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
#include "picongpu/plugins/xrayDiffraction/ReciprocalSpace.hpp"
#include "picongpu/plugins/xrayDiffraction/Writer.hpp"
#include "picongpu/plugins/xrayDiffraction/Implementation.kernel"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/dimensions/DataSpaceOperations.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>
#include <pmacc/memory/MakeUnique.hpp>
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/mpi/reduceMethods/Reduce.hpp>
#include <pmacc/nvidia/functors/Add.hpp>
#include <pmacc/traits/GetNumWorkers.hpp>
#include <pmacc/traits/HasIdentifier.hpp>

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>


namespace picongpu
{
namespace plugins
{
namespace xrayDiffraction
{
namespace detail
{

    template< typename T_Species >
    class Implementation
    {
    public:

        Implementation(
            ReciprocalSpace const & reciprocalSpace,
            std::string const & prefix
        ):
            reciprocalSpace( reciprocalSpace ),
            prefix( prefix )
        {
            isMasterRank = reduce.hasResult(mpi::reduceMethods::Reduce());
            auto totalNumVectors = reciprocalSpace.size.productOfComponents( );
            auto size = DataSpace< DIM1 >( totalNumVectors );
            sumfcoskr = memory::makeUnique< FloatBuffer >( size );
            sumfsinkr = memory::makeUnique< FloatBuffer >( size );
            // allocate one float on GPU and host
            combinedWeighting = memory::makeUnique< FloatBuffer >( DataSpace< DIM1 >( 1 ) );
            // allocate one int on GPU and host
            numMacroparticles = memory::makeUnique< IntBuffer >( DataSpace< DIM1 >( 1 ) );
            globalSumfcoskr.resize( totalNumVectors );
            globalSumfsinkr.resize( totalNumVectors );
            intensity.resize( totalNumVectors );

            if( isMasterRank )
                writer = memory::makeUnique< Writer >( prefix );

        }

        void operator()(
            uint32_t const currentStep,
            MappingDesc const & cellDescription
        )
        {
            computeDiffraction( cellDescription );
            reduceResults();
            if( isMasterRank )
            {
                computeIntensity();
                writer->write(
                    intensity,
                    reciprocalSpace,
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
        void computeDiffraction( MappingDesc const & cellDescription )
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
                KernelXrayDiffraction< numWorkers >{ }
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
} // namespace xrayDiffraction
} // namespace plugins
} // namespace picongpu
