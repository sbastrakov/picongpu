/* Copyright 2020 Sergei Bastrakov
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

#include <picongpu/fields/absorber/Absorber.hpp>
#include <picongpu/fields/antenna/ApplyAntenna.kernel>
#include <picongpu/fields/antenna/Profiles.hpp>

#include <pmacc/math/Vector.hpp>

#include <cstdint>


namespace picongpu
{
namespace fields
{
namespace antenna
{
namespace detail
{

    // Safeguard to not forget to support other antenna types once they are added
    template< typename T_Antenna >
    struct ApplyAntenna
    {
        PMACC_CASSERT_MSG(
            _internal_error_not_supported_antenna_type_used,
            true
        );
    };

    //! Specialization for YMin antenna
    // For now assumes 3d, Yee grid and Yee / YeePML solver
    template<
        typename T_FunctorIncidentE,
        typename T_FunctorIncidentB,
        uint32_t T_gapFromAbsorber
    >
    struct ApplyAntenna<
        YMin<
            T_FunctorIncidentE,
            T_FunctorIncidentB,
            T_gapFromAbsorber
        >
    >
    {
        PMACC_CASSERT_MSG(
            _internal_error_ymin_antenna_only_supported_for_3d,
            simDim == 3
        );

        // This does not work yet because of circular dependencies caused
        /*PMACC_CASSERT_MSG(
            _internal_error_ymin_antenna_only_supported_for_yee_cells,
            std::is_same< Solver::CellType, cellType::Yee >::type::value
        );*/

        // Update E from by using B incident source at t = sourceTimeIteration * dt
        // by delta_t = timeIncrementIteration * dt
        // To be called inside field solver
        void updateE(
            float_X const sourceTimeIteration,
            float_X const timeIncrementIteration,
            MappingDesc const cellDescription
        ) const
        {
            // for now hardcoded for Yee grid like in PIConGPU
            constexpr auto c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
            auto const timeIncrement = timeIncrementIteration * DELTA_T;

            // Ex component uses Bz incident field source
            auto shiftBIncZ = float3_X{ 0.5_X, -0.5_X, 0.0_X };
            // dEx/dt involves c2 * dt * dBz/dy so have to subtract BInc_z
            auto coeffBIncZ = float3_X{ 0.0_X, 0.0_X, -c2 * timeIncrement / cellSize[ 1 ] };

            // Ez component uses Bx incident field source
            auto shiftBIncX = float3_X{ 0.0_X, -0.5_X, 0.5_X };
            // dEz/dt involves -c2 * dt * dBx/dy so have to add BInc_x
            auto coeffBIncX = float3_X{ c2 * timeIncrement / cellSize[ 1 ], 0.0_X, 0.0_X };

            // Create the update functor
            DataConnector & dc = Environment<>::get( ).DataConnector( );
            auto & fieldE = *dc.get< picongpu::FieldE >(
                picongpu::FieldE::getName( ),
                true
            );
            auto dataBox = fieldE.getDeviceDataBox();
            auto const & fieldB = *dc.get< picongpu::FieldB >(
                picongpu::FieldB::getName( ),
                true
            );

            UpdateFunctor<
                T_FunctorIncidentB,
                decltype( dataBox )
            > functor( fieldB.getUnit() );
            functor.field = dataBox;
            functor.inCellShift1 = shiftBIncZ;
            functor.coeff1 = coeffBIncZ;
            functor.inCellShift2 = shiftBIncX;
            functor.coeff2 = coeffBIncX;
            functor.step = sourceTimeIteration;
            // functor.gridIdxShift for now is initialized later
            update(
                functor,
                cellDescription
            );
        }

        // Update B from by using E incident source at t = sourceTimeIteration * dt
        // by delta_t = timeIncrementIteration * dt
        // To be called inside field solver
        void updateB(
            float_X const sourceTimeIteration,
            float_X const timeIncrementIteration,
            MappingDesc const cellDescription
        ) const
        {
            // for now hardcoded for Yee grid like in PIConGPU
            auto const timeIncrement = timeIncrementIteration * DELTA_T;

            // Bx component uses Ez incident field source
            auto shiftEIncZ = float3_X{ 0.0_X, 0.0_X, 0.5_X };
            // dBx/dt involves -dt * dEz/dy so have to add EInc_z
            auto coeffEIncZ = float3_X{ 0.0_X, 0.0_X, timeIncrement / cellSize[ 1 ] };

            // Bz component uses Ex incident field source
            auto shiftEIncX = float3_X{ 0.5_X, 0.0_X, 0.0_X };
            // dBz/dt involves dt * dBx/dy so have to subtract EInc_x
            auto coeffEIncX = float3_X{ -timeIncrement / cellSize[ 1 ], 0.0_X, 0.0_X };

            // Create the update functor
            DataConnector & dc = Environment<>::get( ).DataConnector( );
            auto & fieldB = *dc.get< picongpu::FieldB >(
                picongpu::FieldB::getName( ),
                true
            );
            auto dataBox = fieldB.getDeviceDataBox();
            auto const & fieldE = *dc.get< picongpu::FieldE >(
                picongpu::FieldE::getName( ),
                true
            );

            UpdateFunctor<
                T_FunctorIncidentE,
                decltype( dataBox )
            > functor( fieldE.getUnit() );
            functor.field = dataBox;
            functor.inCellShift1 = shiftEIncZ;
            functor.coeff1 = coeffEIncZ;
            functor.inCellShift2 = shiftEIncX;
            functor.coeff2 = coeffEIncX;
            functor.step = sourceTimeIteration;
            // functor.gridIdxShift for now is initialized later
            update(
                functor,
                cellDescription
            );
        }

    private:

        // time - time at which to take the source.
        // it will be applied to the field at time + dt
        template< typename T_UpdateFunctor >
        void update(
            T_UpdateFunctor functor,
            MappingDesc const cellDescription
        ) const
        {
            // This logic mostly copies laser for now
            auto const step = static_cast< uint32_t >( functor.step );
            bool const antennaNone = false;
            bool const topBoundariesArePeriodic =
               false; /// ( Environment<simDim>::get().GridController().getCommunicationMask( ).isSet( TOP ) );
            const uint32_t numSlides = MovingWindow::getInstance().getSlideCounter(
                step
            );
            bool const boxHasSlided = ( numSlides != 0 );
            bool const disableAntenna =
                antennaNone ||
                topBoundariesArePeriodic ||
                boxHasSlided;
            if( disableAntenna )
                return;

            using Index = pmacc::DataSpace< simDim >;
            using IntVector = pmacc::math::Vector<
                int,
                simDim
            >;
            auto const & subGrid = Environment< simDim >::get().SubGrid();

            // Start and end of the antenna area in the user global coordinates
            // (the coordinate system in which a user functor is expressed,
            // no guards)
            auto beginUserIdx = Index::create( 0 );
            beginUserIdx[ 1 ] += T_gapFromAbsorber;
            auto endUserIdx = subGrid.getGlobalDomain().size;
            for( uint32_t d = 0; d < simDim; d++ )
            {
                beginUserIdx[ d ] += absorber::numCells[ d ][ 0 ];
                endUserIdx[ d ] -= absorber::numCells[ d ][ 1 ];
            }
            // for now hardcode special case for y min
            endUserIdx[ 1 ] = beginUserIdx[ 1 ] + 1;
            /// TODO: remove debug print
            std::cout << "Applying YMin antenna for global y = " << beginUserIdx[ 1 ]
                << " (in the user global coordinate system)\n";

            // Now we express it for the local domain and with grid indices,
            // that include guards. So the indexing to be used in kernels
            auto const localDomain = subGrid.getLocalDomain();
            auto const globalCellOffset = localDomain.offset;
            // Note this conversions so that min and max are happy
            // and we have a DataSpace at the end
            auto const beginLocalUserIdx = Index{
                pmacc::math::max(
                    IntVector{ beginUserIdx - globalCellOffset },
                    IntVector::create( 0 )
                )
            };
            auto const endLocalUserIdx = Index{
                pmacc::math::min(
                    IntVector{ endUserIdx - globalCellOffset },
                    IntVector{ localDomain.size }
                )
            };

            /// TODO: remove debug print
            std::cout << "Local domain offset = " << localDomain.offset
                << ", begin local user idx = " << beginLocalUserIdx << ", end = "
                << endLocalUserIdx << std::endl;

            // Check if we have any active cells in the current domain
            bool areAnyCellsInLocalDomain = true;
            for( uint32_t d = 0; d < simDim; d++ )
                areAnyCellsInLocalDomain = areAnyCellsInLocalDomain &&
                    ( beginLocalUserIdx[ d ] < endLocalUserIdx[ d ] );
            if( !areAnyCellsInLocalDomain )
                return;

            // This setup is copied from laser
            constexpr int laserInitCellsInY = 1;
            using PlaneSizeInSuperCells = typename pmacc::math::CT::AssignIfInRange<
                typename SuperCellSize::vector_type,
                bmpl::integral_c< uint32_t, 1 >, /* y direction */
                bmpl::integral_c< int, laserInitCellsInY >
            >::type;

            auto const superCellSize = SuperCellSize::toRT();
            auto const gridBlocks = ( endLocalUserIdx - beginLocalUserIdx + superCellSize
                - Index::create(1) ) / superCellSize;
            constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                pmacc::math::CT::volume< PlaneSizeInSuperCells >::type::value
            >::value;

            // Shift between internal and user coordinate systems
            pmacc::AreaMapping<
                CORE + BORDER,
                MappingDesc
            > mapper{ cellDescription };
            auto numGuardCells = mapper.getGuardingSuperCells( ) * SuperCellSize::toRT( );

            // Compute corresponding grid indexes and shift to user idx
            auto beginGridIdx = beginLocalUserIdx + numGuardCells;
            auto endGridIdx = endLocalUserIdx + numGuardCells;
            // shift between local grid idx and fractional global user idx:
            // global user idx = local grid idx + idxShift
            functor.gridIdxShift = globalCellOffset - numGuardCells;

            // Kernel call is also analagous to laser
            PMACC_KERNEL(
                ApplyAntennaKernel<
                    numWorkers,
                    PlaneSizeInSuperCells
                >{}
            )(
                gridBlocks,
                numWorkers
            )(
                functor,
                beginGridIdx,
                endGridIdx
            );
        }

    };

    //! Specialization for disabled antenna, just do nothing
    template<>
    struct ApplyAntenna< None >
    {
        void updateE(
            float_X const sourceTimeIteration,
            float_X const timeIncrementIteration,
            MappingDesc const cellDescription
        ) const
        {
        }

        void updateB(
            float_X const sourceTimeIteration,
            float_X const timeIncrementIteration,
            MappingDesc const cellDescription
        ) const
        {
        }
    };

} // namespace detail

    struct ApplyAntenna
    {
        using AntennaSetup = fields::Antenna;
        using AntennaImplementation = detail::ApplyAntenna< AntennaSetup >;

        void updateE(
            float_X const sourceTimeIteration,
            float_X const timeIncrementIteration,
            MappingDesc const cellDescription
        )
        {
            AntennaImplementation{}.updateE(
                sourceTimeIteration,
                timeIncrementIteration,
                cellDescription
            );
        }

        void updateB(
            float_X const sourceTimeIteration,
            float_X const timeIncrementIteration,
            MappingDesc const cellDescription
        )
        {
            AntennaImplementation{}.updateB(
                sourceTimeIteration,
                timeIncrementIteration,
                cellDescription
            );
        }
    };

} // namespace antenna
} // namespace fields
} // namespace picongpu
