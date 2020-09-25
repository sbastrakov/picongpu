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
#include <picongpu/fields/CellType.hpp>
#include <picongpu/fields/Fields.hpp>
#include <picongpu/fields/incidentField/Solver.kernel>
#include <picongpu/fields/incidentField/Profiles.hpp>
#include <picongpu/fields/incidentField/Traits.hpp>

#include <pmacc/math/Vector.hpp>
#include <pmacc/meta/conversion/MakeSeq.hpp>
#include <pmacc/meta/ForEach.hpp>

#include <cstdint>
#include <type_traits>


namespace picongpu
{
namespace fields
{
namespace incidentField
{
namespace detail
{

    struct Parameters
    {
        Parameters(
            MappingDesc const cellDescription
        ):
            cellDescription( cellDescription )
        {
        }

        float_X sourceTimeIteration;

        float_X timeIncrementIteration;

        //! Offset of incident field boundary from min simulation border
        pmacc::DataSpace< simDim > offsetMinBorder;

        //! Offset of incident field boundary from max simulation border
        pmacc::DataSpace< simDim > offsetMaxBorder;

        //! Cell description for kernels
        MappingDesc const cellDescription;

        // +1 or -1
        float_X direction;
    };


    /** Update a field with the given source type
     *
     * @tparam T_Source source type
     */
    template< typename T_Source >
    struct UpdateField;

    /** Safeguard for non-supported source types, instantiating causes a
     *  compile-time error
     *
     * @tparam T_Source source type
     */
    template< typename T_Source >
    struct UpdateField
    {
        PMACC_CASSERT_MSG(
            _internal_error_not_supported_incident_field_source_type_used,
            true
        );
    };
 
    /** Update a field with the incidentField::Source
     *
     * @tparam T_Source source type
     */
    // For now assumes 3d, Yee grid and Yee / YeePML solver
    template<
        typename T_FunctorIncidentE,
        typename T_FunctorIncidentB
    >
    struct UpdateField<
        Source<
            T_FunctorIncidentE,
            T_FunctorIncidentB
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

        // time - time at which to take the source.
        // it will be applied to the field at time + dt
        // UpdatedField += curlCoefficient * Curl(IncidentField)
        template<
            uint32_t T_axis,
            typename T_UpdatedField,
            typename T_IncidentField,
            typename T_IncidentFunctor
        >
        void operator()(
            Parameters const & parameters,
            float_X const curlCoefficient
        )
        {
            // Stop the antenna once the window started moving
            const uint32_t numSlides = MovingWindow::getInstance().getSlideCounter(
                static_cast< uint32_t >( parameters.sourceTimeIteration )
            );
            bool const boxHasSlided = ( numSlides != 0 );
            if( boxHasSlided )
                return;


            /// TODO: figure how to include CellType properly, now there are cyclic dependencies
            auto updatedFieldPositions = traits::FieldPosition<
                cellType::Yee,
                T_UpdatedField
            >{}();
            auto incidentFieldPositions = traits::FieldPosition<
                cellType::Yee,
                T_IncidentField
            >{}();

            /* Whether the field values we are updating are in the TF or SF region:
             * TF when x component is not on the x cell border,
             * SF when x component is on the x cell border
             */
            bool isUpdatedFieldTF = ( updatedFieldPositions[0][0] != 0.0_X );

            using Index = pmacc::DataSpace< simDim >;
            using IntVector = pmacc::math::Vector<
                int,
                simDim
            >;
            auto const & subGrid = Environment< simDim >::get().SubGrid();

            // Start and end of the antenna area in the user global coordinates
            // (the coordinate system in which a user functor is expressed,
            // no guards)
            auto beginUserIdx = parameters.offsetMinBorder;
            // TF in positive direction needs to have the begin adjusted by one to get
            // the first index what's inside the TF region
            if( isUpdatedFieldTF && parameters.direction > 0 )
                beginUserIdx += pmacc::DataSpace< simDim >::create( 1 );
            auto endUserIdx = subGrid.getGlobalDomain().size - parameters.offsetMaxBorder;
            if( parameters.direction > 0 )
                endUserIdx[ T_axis ] = beginUserIdx[ T_axis ] + 1;
            else
                beginUserIdx[ T_axis ] = endUserIdx[ T_axis ] - 1;

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

            // Check if we have any active cells in the current domain
            bool areAnyCellsInLocalDomain = true;
            for( uint32_t d = 0; d < simDim; d++ )
                areAnyCellsInLocalDomain = areAnyCellsInLocalDomain &&
                ( beginLocalUserIdx[ d ] < endLocalUserIdx[ d ] );
            if( !areAnyCellsInLocalDomain )
                return;

            // This setup is copied from laser
            using PlaneSizeInSuperCells = typename pmacc::math::CT::AssignIfInRange<
                typename SuperCellSize::vector_type,
                bmpl::integral_c< uint32_t, T_axis >,
                bmpl::integral_c< int, 1 >
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
            > mapper{ parameters.cellDescription };
            auto numGuardCells = mapper.getGuardingSuperCells( ) * SuperCellSize::toRT( );

            // Compute corresponding grid indexes and shift to user idx
            auto beginGridIdx = beginLocalUserIdx + numGuardCells;
            auto endGridIdx = endLocalUserIdx + numGuardCells;

            // dir0 is boundary normal axis, dir1 and dir2 are two other axes
            constexpr auto dir0 = T_axis;
            constexpr auto dir1 = ( dir0 + 1 ) % 3;
            constexpr auto dir2 = ( dir0 + 2 ) % 3;

            auto onAxisShift = float3_X::create( 0.0_X );
            float_X const incidentFieldSign = 0.0_X;
            if( isUpdatedFieldTF )
            {
                /* When updating the total field, the incident field needs to be
                 * added to SF terms for external neighbors which are half-cell shifted to the outside
                 */
                onAxisShift[ dir0 ] = -0.5_X * parameters.direction;
            }
            else
            {
                /* When updating the scattered field, the incident field needs to be
                 * subtracted from TF terms for internal neighbors which are half-cell shifted to the inside
                 */
                onAxisShift[ dir0 ] = 0.5_X * parameters.direction;
            }

            // The shifts are two cross-combinations of dir1, dir2 components

            auto shiftInc1 = incidentFieldPositions[ dir1 ] + onAxisShift;
            auto shiftInc2 = incidentFieldPositions[ dir2 ] + onAxisShift;

            float_X const coeffBase = incidentFieldSign * curlCoefficient /
                cellSize[ T_axis ];
            auto coeffInc1 = float3_X::create( 0.0_X );
            coeffInc1[ dir1 ] = coeffBase;
            auto coeffInc2 = float3_X::create( 0.0_X );
            coeffInc2[ dir2 ] = -coeffBase;

            // For the propagation direction the shift is half cell outwards
            // Create the update functor
            DataConnector & dc = Environment<>::get( ).DataConnector( );
            auto & updatedField = *dc.get< T_UpdatedField >(
                T_UpdatedField::getName( ),
                true
            );
            auto dataBox = updatedField.getDeviceDataBox();
            auto const & incidentField = *dc.get< T_IncidentField >(
                T_IncidentField::getName( ),
                true
            );

            UpdateFunctor<
                T_IncidentFunctor,
                decltype( dataBox )
            > functor( incidentField.getUnit() );
            functor.field = dataBox;
            // take care of 2d
            for( uint32_t d = 0; d < simDim; d++ )
            {
                functor.inCellShift1[ d ] = shiftInc1[ d ];
                functor.inCellShift2[ d ] = shiftInc2[ d ];
            }
            functor.coeff1 = coeffInc1;
            functor.coeff2 = coeffInc2;
            functor.step = parameters.sourceTimeIteration;
            // shift between local grid idx and fractional global user idx:
            // global user idx = local grid idx + idxShift
            functor.gridIdxShift = globalCellOffset - numGuardCells;

            std::cout << "begin local user idx = " << beginLocalUserIdx << ", end = "
                << endLocalUserIdx << std::endl;
            std::cout << "inCellShift1 = " << functor.inCellShift1 << ", coeff1 = " <<
                functor.coeff1 << "\n";
            std::cout << "inCellShift2 = " << functor.inCellShift2 << ", coeff1 = " <<
                functor.coeff2 << "\n";

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


    template< typename T_Source >
    struct UpdateE
    {
        PMACC_CASSERT_MSG(
            _internal_error_not_supported_incident_field_source_type_used,
            true
        );
    };

    //! Specialization for None source does nothing
    template< >
    struct UpdateE< None >
    {
        template< uint32_t T_axis >
        void operator()( Parameters const & parameters )
        {
        }
    };

    template<
        typename T_FunctorIncidentE,
        typename T_FunctorIncidentB
    >
    struct UpdateE<
        Source<
            T_FunctorIncidentE,
            T_FunctorIncidentB
        >
    >
    {
        // Update E from by using B incident source at t = sourceTimeIteration * dt
        // by delta_t = timeIncrementIteration * dt
        // To be called inside field solver
        template< uint32_t T_axis >
        void operator()( Parameters const & parameters )
        {
            /* E field values are located on the internal side of the Huygence
             * surface, so they are TF, neighbor internal B values are TF and
             * the neighbor external B values are SF.
             *
             * Therefore, the finite difference expression for such values of E
             * should result in a TF. The field solver update used TF for
             * neighbor internal B and SF for neighbor external B. Thus, here we
             * need to add incident B values to terms for neighbor external B in
             * the finite difference derivatives along T_axis.
             *
             * For the max border such B neighbors are the same grid index,
             * while for the min border the B neighbors correspond to the
             * Huygence surface offset, but the E grid indices are one more.
             */
            using UpdatedField = picongpu::FieldE;
            using IncidentField = picongpu::FieldB;
            constexpr auto c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
            auto const timeIncrement = parameters.timeIncrementIteration *
                DELTA_T;
            // E(t + timeIncrement) = E(t) + timeIncrement * c2 * curl(B)
            auto const curlCoefficient = timeIncrement * c2;
            UpdateField<
                Source<
                    T_FunctorIncidentE,
                    T_FunctorIncidentB
                >
            > updateField;
            updateField.template operator()<
                T_axis,
                UpdatedField,
                IncidentField,
                T_FunctorIncidentB
            >(
                parameters,
                curlCoefficient
            );
        }
    };

    template< typename T_Source >
    struct UpdateB
    {
        PMACC_CASSERT_MSG(
            _internal_error_not_supported_incident_field_source_type_used,
            true
        );
    };

    //! Specialization for None source does nothing
    template< >
    struct UpdateB< None >
    {
        template< uint32_t T_axis >
        void operator()( Parameters const & parameters )
        {
        }
    };

    template<
        typename T_FunctorIncidentE,
        typename T_FunctorIncidentB
    >
    struct UpdateB<
        Source<
            T_FunctorIncidentE,
            T_FunctorIncidentB
        >
    >
    {
        // Update B from by using E incident source at t = sourceTimeIteration * dt
        // by delta_t = timeIncrementIteration * dt
        // To be called inside field solver
        template< uint32_t T_axis >
        void operator()( Parameters const & parameters )
        {
            using UpdatedField = picongpu::FieldB;
            using IncidentField = picongpu::FieldE;
            auto const timeIncrement = parameters.timeIncrementIteration *
                DELTA_T;
            // B(t + timeIncrement) = B(t) - timeIncrement * curl(E)
            auto const curlCoefficient = -timeIncrement;
            UpdateField<
                Source<
                    T_FunctorIncidentE,
                    T_FunctorIncidentB
                >
            > updateField;
            updateField.template operator()<
                T_axis,
                UpdatedField,
                IncidentField,
                T_FunctorIncidentE
            >(
                parameters,
                curlCoefficient
            );
        }
    };

} // namespace detail

    /** Solver for incident fields
     *
     * Incident fields are generated by the virtual Huygens surface.
     * The Huygens surface is geometrically parallel to the interface between the
     * absorber area and internal (non-absorbing) simulation area. Thus, it is
     * composed of six axis-aligned plane segments in 3D and four axis-aligned
     * line segments in 2D.
     * The position of the Huygens surface is controlled by offset from the
     * absorber-intenral area interface. This offset is a sum of the
     * user-specified gap and 0.75 of a cell into the internal area from each
     * side. (The choice of 0.75 is arbitrary, it is only required that the
     * surface is not aligned with full or half cells). The implementation
     * provides correct treatment of any gap >= 0 as long as the surface has
     * at least one cell in the internal volume.
     */
    struct Solver
    {

        /** Z sources to be used by implementation, overrides the names from
         *  incidentField.param on purpose to handle 2d and 3d cases uniformly
         */
        using ZMin = GetZMin::type;
        using ZMax = GetZMax::type;

        //! Incident fields for all borders
        using AllIncidentFields = pmacc::MakeSeq_t<
            XMin,
            XMax,
            YMin,
            YMax,
            ZMin,
            ZMax
        >;

        Solver( MappingDesc const cellDescription ):
            cellDescription( cellDescription )
        {
            // Compute offsets from local domain borders, without guards
            absorber::Thickness thickness = absorber::getGlobalThickness();
            for( uint32_t axis = 0u; axis < simDim; axis++ )
            {
                offsetMinBorder[ axis ] = thickness( axis, 0 ) +
                    GAP_FROM_ABSORBER[ axis ][ 0 ];
                offsetMaxBorder[ axis ] = thickness( axis, 1 ) +
                    GAP_FROM_ABSORBER[ axis ][ 1 ];
            }
        }

        void updateE( float_X const sourceTimeIteration )
        {
            auto parameters = detail::Parameters{
                cellDescription
            };
            parameters.sourceTimeIteration = sourceTimeIteration;
            parameters.timeIncrementIteration = 1.0_X;
            parameters.offsetMinBorder = offsetMinBorder;
            parameters.offsetMaxBorder = offsetMaxBorder;

            parameters.direction = 1.0_X;
            std::cout << "\n\n\n\nUpdateE Xmin\n";
            detail::UpdateE< XMin >{}.operator()< 0 >( parameters );
            parameters.direction = -1.0_X;
            std::cout << "\n\nUpdateE Xmax\n";
            detail::UpdateE< XMax >{}.operator()< 0 >( parameters );

            parameters.direction = 1.0_X;
            std::cout << "\n\nUpdateE Ymin\n";
            detail::UpdateE< YMin >{}.operator()< 1 >( parameters );
            parameters.direction = -1.0_X;
            std::cout << "\nUpdateE Ymax\n";
            detail::UpdateE< YMax >{}.operator()< 1 >( parameters );

            parameters.direction = 1.0_X;
            std::cout << "\n\nUpdateE Zmin\n";
            detail::UpdateE< ZMin >{}.operator()< 2 >( parameters );
            parameters.direction = -1.0_X;
            std::cout << "\nUpdateE Zmax\n";
            detail::UpdateE< ZMax >{}.operator()< 2 >( parameters );
        }

        void updateBHalf( float_X const sourceTimeIteration )
        {
            auto parameters = detail::Parameters{
                cellDescription
            };
            parameters.sourceTimeIteration = sourceTimeIteration;
            parameters.timeIncrementIteration = 0.5_X;
            parameters.offsetMinBorder = offsetMinBorder;
            parameters.offsetMaxBorder = offsetMaxBorder;

            parameters.direction = 1.0_X;
            std::cout << "\n\n\n\nUpdateB Xmin\n";
            detail::UpdateB< XMin >{}.operator()< 0 >( parameters );
            parameters.direction = -1.0_X;
            std::cout << "\n\nUpdateB Xmax\n";
            detail::UpdateB< XMax >{}.operator()< 0 >( parameters );

            parameters.direction = 1.0_X;
            std::cout << "\n\nUpdateB Ymin\n";
            detail::UpdateB< YMin >{}.operator()< 1 >( parameters );
            parameters.direction = -1.0_X;
            std::cout << "\nUpdateB Ymax\n";
            detail::UpdateB< YMax >{}.operator()< 1 >( parameters );

            parameters.direction = 1.0_X;
            std::cout << "\n\nUpdateB Zmin\n";
            detail::UpdateB< ZMin >{}.operator()< 2 >( parameters );
            parameters.direction = -1.0_X;
            std::cout << "\nUpdateB Zmax\n";
            detail::UpdateB< ZMax >{}.operator()< 2 >( parameters );
        }

    private:

        /** Offset of the Hyugence surface from min borders of the global domain
         *
         * The offset is from the start of CORE + BORDER area.
         * Counted in full cells, the surface is additionally offset by 0.75 cells
         */
        pmacc::DataSpace< simDim > offsetMinBorder;

        /** Offset of the Hyugence surface from max borders of the global domain
         *
         * The offset is from the end of CORE + BORDER area.
         * Counted in full cells, the surface is additionally offset by 0.75 cells
         */
        pmacc::DataSpace< simDim > offsetMaxBorder;

        //! Cell description for kernels
        MappingDesc const cellDescription;

    };

} // namespace incidentField
} // namespace fields
} // namespace picongpu
