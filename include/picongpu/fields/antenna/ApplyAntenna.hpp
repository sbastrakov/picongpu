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
///#include <picongpu/fields/MaxwellSolver/YeePML/Parameters.hpp>

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
        template< typename T_FieldJ >
        void operator()(
            T_FieldJ & fieldJ,
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
    template< typename T_Functor >
    struct ApplyAntenna< YMin< T_Functor > >
    {
        template< typename T_FieldJ >
        void operator()(
            T_FieldJ & fieldJ,
            uint32_t const step,
            MappingDesc const cellDescription
        ) const
        {
            // This logic mostly copies laser for now
            bool const antennaNone = false;
            bool const topBoundariesArePeriodic =
                ( Environment<simDim>::get().GridController().getCommunicationMask( ).isSet( TOP ) );
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
            // shift between local grid idx and global user idx:
            // global user idx = local grid idx + idxShift
            auto const idxShift = globalCellOffset - numGuardCells;

            auto functor = typename YMin< T_Functor >::Functor{ fieldJ.getUnit() };

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
                fieldJ.getDeviceDataBox( ),
                step,
                beginGridIdx,
                endGridIdx,
                idxShift
            );
        }

    private:

        //using Thickness = maxwellSolver::yeePML::Thickness;

        ///** Get absorber thickness for the local domain at the current time step.
        // *
        // * It depends on the current step because of the moving window.
        // * Currently the function is almost a copy of yeePML::Solver::getLocalThickness()
        // * Perhaps this logic should be moved to general absorber and reused here and from PML
        // */
        //Thickness getLocalThickness( uint32_t const currentStep ) const
        //{
        //    /* The logic of the following checks is the same as in
        //    * absorber::ExponentialDamping::run( ), to disable the absorber
        //    * at a border we set the corresponding thickness to 0.
        //    */
        //    auto & movingWindow = MovingWindow::getInstance( );
        //    auto const numSlides = movingWindow.getSlideCounter( currentStep );
        //    auto const numExchanges = NumberOfExchanges< simDim >::value;
        //    auto const communicationMask = Environment< simDim >::get( ).GridController( ).getCommunicationMask( );
        //    // First set to global absorber thickness
        //    Thickness thickness;
        //    for( uint32_t axis = 0u; axis < simDim; axis++ )
        //        for( auto direction = 0; direction < 2; direction++ )
        //            thickness( axis, direction ) = absorber::numCells[ axis ][ direction ];
        //    for( uint32_t exchange = 1u; exchange < numExchanges; ++exchange )
        //    {
        //        /* Here we are only interested in the positive and negative
        //        * directions for x, y, z axes and not the "diagonal" ones.
        //        * So skip other directions except left, right, top, bottom,
        //        * back, front
        //        */
        //        if( FRONT % exchange != 0 )
        //            continue;

        //        // Transform exchange into a pair of axis and direction
        //        uint32_t axis = 0;
        //        if( exchange >= BOTTOM && exchange <= TOP )
        //            axis = 1;
        //        if( exchange >= BACK )
        //            axis = 2;
        //        uint32_t direction = exchange % 2;

        //        // No PML at the borders between two local domains
        //        bool hasNeighbour = communicationMask.isSet( exchange );
        //        if( hasNeighbour )
        //            thickness( axis, direction ) = 0;

        //        // Disable PML at the far side of the moving window
        //        if( movingWindow.isSlidingWindowActive( currentStep ) && exchange == BOTTOM )
        //            thickness( axis, direction ) = 0;
        //    }
        //    return thickness;
        //}

    };

    //! Specialization for disabled antenna, just do nothing
    template<>
    struct ApplyAntenna< None >
    {
        template< typename T_FieldJ >
        void operator()(
            T_FieldJ & fieldJ,
            uint32_t const step,
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
