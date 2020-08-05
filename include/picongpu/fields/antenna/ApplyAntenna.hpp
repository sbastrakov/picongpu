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
#include <picongpu/fields/CellType.hpp>
#include <picongpu/traits/FieldPosition.hpp>

#include <pmacc/math/Vector.hpp>

#include <cstdint>


namespace picongpu
{
namespace fields
{
namespace antenna
{

    template< typename T_Antenna >
    struct ApplyAntenna
    {
        template< typename T_Field >
        void operator()(
            T_Field & field,
            uint32_t const step,
            MappingDesc const cellDescription
        ) const
        {
            PMACC_CASSERT_MSG(
                _internal_error_not_implemented_antenna_used,
                true
            );
        }
    };

    //! Specialization for YMin antenna
    template<
        typename T_FunctorE,
        typename T_FunctorB,
        uint32_t T_gapFromAbsorber
    >
    struct ApplyAntenna<
        YMin<
            T_FunctorE,
            T_FunctorB,
            T_gapFromAbsorber
        >
    >
    {
        // time - time at which to take the source.
        // it will be applied to the field at time + dt
        template< typename T_Field >
        void operator()(
            T_Field & field,
            float_X const step,
            MappingDesc const cellDescription
        ) const
        {
            // This logic mostly copies laser for now
            bool const antennaNone = false;
            bool const topBoundariesArePeriodic =
                ( Environment<simDim>::get().GridController().getCommunicationMask( ).isSet( TOP ) );
            const uint32_t numSlides = MovingWindow::getInstance().getSlideCounter(
                static_cast< uint32_t >( step )
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

            // shift inside the cell, in units of cells
            //// ISSUE: this is vector of vectors, see notes
            using FieldPosition = typename traits::FieldPosition<
                fields::CellType,
                T_Field
            >;
            auto const shiftInCell = FieldPosition{}();

            // Compute corresponding grid indexes and shift to user idx
            auto beginGridIdx = beginLocalUserIdx + numGuardCells;
            auto endGridIdx = endLocalUserIdx + numGuardCells;
            // shift between local grid idx and fractional global user idx:
            // global user idx = local grid idx + idxShift
            auto const idxShift = pmacc::algorithms::precisionCast::precisionCast< float_X >(
                globalCellOffset - numGuardCells ) + shiftInCell;

            // TODO: generalize
            auto functor = T_FunctorE{ field.getUnit() };

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
                field.getDeviceDataBox( ),
                step,
                beginGridIdx,
                endGridIdx,
                idxShift
            );
        }

    };

    //! Specialization for disabled antenna, just do nothing
    template<>
    struct ApplyAntenna< None >
    {
        template< typename T_Field >
        void operator()(
            T_Field & field,
            float_X const time,
            MappingDesc const cellDescription
        ) const
        {
            /// TODO: remove debug print
            std::cout << "Applying none antenna\n";
        }
    };

} // namespace antenna
} // namespace fields
} // namespace picongpu
