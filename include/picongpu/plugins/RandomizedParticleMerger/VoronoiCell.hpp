/* Copyright 2017-2020 Heiko Burau, Xeinia Bastrakova, Sergei Bastrakov
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

#include "picongpu/algorithms/KinEnergy.hpp"

#include <pmacc/types.hpp>

#include <cstdint>


namespace picongpu
{
namespace plugins
{
namespace randomizedParticleMerger
{

    //! Status of a Voronoi cell
    enum struct VoronoiStatus : uint8_t
    {
        /* !< a Voronoi cell is collecting particles (first state) */
        collecting,
        /* !< the Voronoi cell is splitting thus all its particles have
         * to move to one of two sub-Voronoi cells
         */
        splitting,
        /* !< the cell needs to be destroyed. Before this can happen
         * all its particles need to clear their voronoiCellId attribute.
         */
        abort,
        /* !< the Voronoi cell is ready for merging. After merging it is destroyed. */
        readyForMerging,
    };


    /** Stage of a Voronoi cell
     *
     * The spliiting process is two-fold: at first, the splitting is done regarding
     * only the spread in position and then by looking at the spread of momentum.
     */
    enum struct VoronoiSplittingStage : bool
    {
        /* !< the spatial distribution is splitted */
        position,
        /* !< the momentum distribution is splitted */
        momentum
    };

    //! Voronoi cell representation
    struct VoronoiCell
    {
        VoronoiStatus status;
        VoronoiSplittingStage splittingStage;
        uint32_t numMacroParticles;
        float_X numRealParticles;

        float3_X meanMomentumValue;
        float3_X meanPositionValue;
        float3_X meanMomentumSquaredValue;
        float3_X meanPositionSquaredValue;

        uint8_t splittingComponent;
        int32_t lowerCellId;
        int32_t higherCellId;
        int firstParticleFlag; // TODO: int32_t ? beware of atomics; keep int if complicated
        float_X expectedNumMacroParticles;
        uint32_t parentNumMacroParticles;
        float_X parentExpectedNumMacroParticles;

        HDINLINE
        VoronoiCell(
            VoronoiSplittingStage splittingStage = VoronoiSplittingStage::position,
            float_X parentNumMacroParticles = 0.0_X,
            float_X parentExpectedNumMacroParticles = -1.0_X
        ) :
            status( VoronoiStatus::collecting ),
            splittingStage( splittingStage ),
            numMacroParticles( 0u ),
            numRealParticles( float_X( 0.0_X ) ),
            meanMomentumValue( float3_X::create( 0.0_X ) ),
            meanPositionValue( float3_X::create( 0.0_X ) ),
            meanMomentumSquaredValue( float3_X::create( 0.0_X ) ),
            meanPositionSquaredValue( float3_X::create( 0.0_X ) ),
            firstParticleFlag( 0 ),
            expectedNumMacroParticles ( 0 ),
            parentNumMacroParticles ( parentNumMacroParticles ),
            parentExpectedNumMacroParticles ( parentExpectedNumMacroParticles )

        {}

        /** status setter */
        HDINLINE
        void setToAbort()
        {
            status = VoronoiStatus::abort;
        }


        /** Mark the cell for splitting
         *
         * @param splittingComponent index of position or momentum component
         *                           to use for splitting
         * @param lowerCellId cell index of a new "lower" subcell
         * @param higherCellId cell index of a new "upper" subcell
         */
        HDINLINE
        void setToSplitting(
            const uint8_t splittingComponent,
            const int32_t lowerCellId,
            const int32_t higherCellId)
        {
            status = VoronoiStatus::splitting;
            this->splittingComponent = splittingComponent;
            this->lowerCellId = lowerCellId;
            this->higherCellId = higherCellId;
        }


        /** status setter */
        HDINLINE
        void setToReadyForMerging()
        {
            this->status = VoronoiStatus::readyForMerging;
        }

        /** check if the current thread is associated to the first particle */
        template< typename T_Acc >
        DINLINE
        bool isFirstParticle(T_Acc const & acc)
        {
            return atomicExch( &this->firstParticleFlag, 1 ) == 0;
        }


        /** add a particle to this Voronoi cell */
        template< typename T_Acc >
        DINLINE
        void addParticle(
            T_Acc const & acc,
            const floatD_X position,
            const float3_X momentum,
            const float_X weighting
        )
        {
            // TODO: move to modern cupla
            atomicAdd( &this->numMacroParticles, (uint32_t)1, ::alpaka::hierarchy::Threads{} );
            atomicAdd( &this->numRealParticles, weighting, ::alpaka::hierarchy::Threads{} );

            const floatD_X position2 = position * position;

            for( int i = 0; i < simDim; i++ )
            {
                atomicAdd( &this->meanPositionValue[i], weighting * position[i], ::alpaka::hierarchy::Threads{} );
                atomicAdd( &this->meanPositionSquaredValue[i], weighting * position2[i], ::alpaka::hierarchy::Threads{} );
            }

            const float3_X momentum2 = momentum * momentum;

            for( int i = 0; i < DIM3; i++ )
            {
                atomicAdd( &this->meanMomentumValue[i], weighting * momentum[i], ::alpaka::hierarchy::Threads{} );
                atomicAdd( &this->meanMomentumSquaredValue[i], weighting * momentum2[i], ::alpaka::hierarchy::Threads{} );
            }
        }

        /** Counting parameters that are necessary before processing vornoi cell:
         *  mean values and expected number of macro particles
         *
         * @param minMacroParticlesToDivide min number of macroparticles in a cell
         *                                  such that the cell is always subdivided
         * @param ratioKeptParticles ratio of particles that are kept on average
         */
        HDINLINE
        void finalizePrecalculationValues(
            const uint32_t minMacroParticlesToDivide,
            const float_X ratioKeptParticles
        )
        {
            finalizeMeanValues();
            finalizeExpectedNumberParticles( minMacroParticlesToDivide, ratioKeptParticles );

        }

        //! Finalize calculation of mean values
        HDINLINE
        void finalizeMeanValues()
        {
            meanMomentumValue /= numRealParticles;
            meanPositionValue /= numRealParticles;
            meanMomentumSquaredValue /= numRealParticles;
            meanPositionSquaredValue /= numRealParticles;
        }

        /** Count expected number of particles in the cell
         *
         * @param minMacroParticlesToDivide min number of macroparticles in a cell
         *                                  such that the cell is always subdivided
         * @param ratioKeptParticles ratio of particles that are kept on average
         */
        HDINLINE
        void finalizeExpectedNumberParticles(
            uint32_t const minMacroParticlesToDivide,
            float_X const ratioKeptParticles
        )
        {
            // Special case for the original voronoi cells
            if ( parentExpectedNumMacroParticles < 0 )
            {
                expectedNumMacroParticles = numMacroParticles * ratioKeptParticles;
                return;
            }

            // Algorithm stop conditions for 1 and 2 macroparticles
            if ( numMacroParticles == 1 )
                expectedNumMacroParticles = 1;
            if ( numMacroParticles == 2 && parentNumMacroParticles == 3 )
                expectedNumMacroParticles = 2;

            // Normal subdivision step
            if ( parentNumMacroParticles > minMacroParticlesToDivide )
            {
                expectedNumMacroParticles = numMacroParticles * ratioKeptParticles;
            }
            else
            {
                // TODO: this equation works, but there may be a better one (see the notebook)
                float_X undividedCellCoeff = (
                    parentExpectedNumMacroParticles + parentNumMacroParticles ) / 2.0_X;
                float_X currentExpectedNumMacroParticles =
                    numMacroParticles * undividedCellCoeff / parentNumMacroParticles;
                expectedNumMacroParticles = currentExpectedNumMacroParticles;
            }
        }

        /** determine in which of the two sub-Voronoi cells a particle falls */
        HDINLINE
        int32_t getSubVoronoiCell(
            const floatD_X position,
            const float3_X momentum
        ) const
        {
            const float_X valParticle =
                splittingStage == VoronoiSplittingStage::position ?
                position[splittingComponent] :
                momentum[splittingComponent];
            const float_X meanVoronoi =
                splittingStage == VoronoiSplittingStage::position ?
                meanPositionValue[splittingComponent] :
                meanMomentumValue[splittingComponent];
            return
                valParticle < meanVoronoi ?
                lowerCellId :
                higherCellId;
        }

        /** auxillary function for getting the mean squared deviation in position or momentum */
        HDINLINE
        float_X getMaxValueSpread2(
            uint8_t& component,
            const uint8_t dimension
        ) const
        {
            const float3_X meanValue2 =
                  splittingStage == VoronoiSplittingStage::position ?
                  meanPositionValue * meanPositionValue :
                  meanMomentumValue * meanMomentumValue;

            const float3_X valueSpread2 =
                  splittingStage == VoronoiSplittingStage::position ?
                  meanPositionSquaredValue - meanValue2 :
                  meanMomentumSquaredValue - meanValue2;

            /* find component of most spread in position */
            component = 0;
            float_X maxValueSpread2 = valueSpread2[0];
            for( uint8_t i = 1; i < dimension; i++ )
            {
                if( valueSpread2[i] > maxValueSpread2 )
                {
                    maxValueSpread2 = valueSpread2[i];
                    component = i;
                }
            }

            return maxValueSpread2;
        }


        /** calculate the maxmimum squared spread in position
         *
         * @param component index of position component of maxmimum spread
         * @return maxmimum squared spread in position
         */
        HDINLINE
        float_X getMaxPositionSpread2( uint8_t& component ) const
        {
            return getMaxValueSpread2( component, simDim );
        }


        /** calculate the maxmimum squared spread in momentum
         *
         * @param component index of momentum component of maxmimum spread
         * @return maxmimum squared spread in momentum
         */
        HDINLINE
        float_X getMaxMomentumSpread2( uint8_t& component ) const
        {
            return getMaxValueSpread2( component, DIM3 );
        }

        /** invesing splitting stage */
        HDINLINE
        void invertSplittingStage( )
        {
            if( splittingStage == VoronoiSplittingStage::position )
                splittingStage = VoronoiSplittingStage::momentum;
            else
                splittingStage = VoronoiSplittingStage::position;
        }
    };

} // namespace randomizedParticleMerger
} // namespace plugins
} // namespace picongpu
