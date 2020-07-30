/* Copyright 2013-2020 Axel Huebl, Heiko Burau, Rene Widera, Remi Lehe
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

#include "picongpu/fields/MaxwellSolver/Lehe/Derivative.def"
#include "picongpu/fields/differentiation/Derivative.hpp"
#include "picongpu/fields/differentiation/ForwardDerivative.hpp"
#include "picongpu/fields/differentiation/Traits.hpp"

#include <pmacc/algorithms/math/defines/pi.hpp>
#include <pmacc/types.hpp>
#include <pmacc/math/Vector.hpp>

#include "picongpu/traits/GetMargin.hpp"

#include <pmacc/meta/accessors/Identity.hpp>


namespace picongpu
{
namespace fields
{
namespace maxwellSolver
{
namespace lehe
{

    template<
        uint32_t T_cherenkovFreeDirection,
        uint32_t T_direction
    >
    struct DerivativeFunctor;

    // Differentiate at the T_cherenkovFreeDirection
    template<
        uint32_t T_direction
    >
    struct DerivativeFunctor<
        T_direction,
        T_direction
    >
    {
        float_X mySin;

        HINLINE DerivativeFunctor( )
        {
            mySin = float_X(
                math::sin(
                    pmacc::math::Pi< float_64 >::halfValue *
                    float_64( SPEED_OF_LIGHT ) * float_64( DELTA_T ) /
                    float_64( cellSize[ T_direction ] )
                )
            );
        }

        template< typename T_DataBox >
        HDINLINE typename T_DataBox::ValueType operator( )(
            T_DataBox const & data ) const
        {
            constexpr uint32_t dir0 = T_direction;
            constexpr uint32_t dir1 = (dir0 + 1) % 3;
            constexpr uint32_t dir2 = (dir0 + 2) % 3;

            constexpr float_X dt2 = DELTA_T * DELTA_T;
            constexpr float_X c2 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;

            // cellSize is not constexpr currently, so make an own constexpr array
            constexpr float_X size[3] = { CELL_WIDTH, CELL_HEIGHT, CELL_DEPTH };
            constexpr float_X betaDir1 = 0.125_X * size[ dir0 ] * size[ dir0 ]
                / (size[ dir1 ] * size[ dir1 ]);
            constexpr float_X betaDir2 = 0.125_X * size[ dir0 ] * size[ dir0 ]
                / (size[ dir1 ] * size[ dir1 ]);

            float_X const delta = 0.25_X *
                ( 1.0_X - cellSize[ dir0 ] * cellSize[ dir0 ] / ( c2 * dt2 ) * mySin * mySin );

            float_X const alpha = 1.0_X - 2.0_X * betaDir1 - 2.0_X * betaDir2
                - 3.0_X * delta;

            pmacc::DataSpace< simDim > const currentIndex;
            pmacc::DataSpace< simDim > secondUpperIndexDir0, lowerIndexDir0;
            secondUpperIndexDir0[ dir0 ] = 2;
            lowerIndexDir0[ dir0 ] = -1;
            pmacc::DataSpace< simDim > upperNeighborDir1, lowerNeighborDir1;
            upperNeighborDir1[ dir1 ] = 1;
            lowerNeighborDir1[ dir1 ] = -1;
            pmacc::DataSpace< simDim > upperNeighborDir2, lowerNeighborDir2;
            upperNeighborDir2[ dir2 ] = 1;
            lowerNeighborDir2[ dir2 ] = -1;
            auto forwardDerivative = differentiation::makeDerivativeFunctor<
                differentiation::Forward,
                T_direction
            >();
            return
                alpha * forwardDerivative( data ) +
                betaDir1 * forwardDerivative( data.shift( upperNeighborDir1 ) ) +
                betaDir1 * forwardDerivative( data.shift( lowerNeighborDir1 ) ) +
                betaDir2 * forwardDerivative( data.shift( upperNeighborDir2 ) ) +
                betaDir2 * forwardDerivative( data.shift( lowerNeighborDir2 ) ) +
                delta * ( data( secondUpperIndexDir0 ) - data( lowerIndexDir0 ) ) /
                cellSize[ T_direction ];
        }
    };

    // Differentiate at another direction
    template<
        uint32_t T_cherenkovFreeDirection,
        uint32_t T_direction
    >
    struct DerivativeFunctor
    {
        PMACC_CASSERT_MSG(
            _internal_error_wrong_lehe_derivative_functor_specialization,
            T_cherenkovFreeDirection != T_direction
        );

        template< typename T_DataBox >
        HDINLINE typename T_DataBox::ValueType operator( )(
            T_DataBox const & data ) const
        {
            /* In this case equation (6) greatly simplifies with
             * delta = 0 and one of the two beta coefficients is also 0.
             * So we only need the simple forward derivative operator and take
             * the weighted average of its values at the current point and the
             * two neighbors in T_cherenkovFreeDirection
             */
            constexpr float_X beta = 0.125_X;
            constexpr float_X alpha = 1.0_X - 2.0_X * beta;
            auto forwardDerivative = differentiation::makeDerivativeFunctor<
                differentiation::Forward,
                T_direction
            >();
            pmacc::DataSpace< simDim > const currentIndex;
            pmacc::DataSpace< simDim > upperNeighbor;
            upperNeighbor[ T_cherenkovFreeDirection ] = 1;
            pmacc::DataSpace< simDim > lowerNeighbor;
            lowerNeighbor[ T_cherenkovFreeDirection ] = -1;
            return alpha * forwardDerivative( data ) +
                beta * forwardDerivative( data.shift( upperNeighbor ) ) +
                beta * forwardDerivative( data.shift( lowerNeighbor ) );
        }
    };

} // namespace lehe
} // namespace maxwellSolver


namespace differentiation
{
namespace traits
{

    /** Functor type trait specialization for forward derivative
     *
     * @tparam T_direction direction to take derivative in, 0 = x, 1 = y, 2 = z
     */
    template<
        uint32_t T_cherenkovFreeDirection,
        uint32_t T_direction
    >
    struct DerivativeFunctor<
        maxwellSolver::lehe::Derivative< T_cherenkovFreeDirection >,
        T_direction
    > : pmacc::meta::accessors::Identity<
        maxwellSolver::lehe::DerivativeFunctor<
            T_cherenkovFreeDirection,
            T_direction
        >
    >
    {
    };

} // namespace traits
} // namespace differentiation
} // namespace fields

namespace traits
{

    //! Get margin of the forward derivative
    template< uint32_t T_cherenkovFreeDirection >
    struct GetMargin<
        fields::maxwellSolver::lehe::Derivative< T_cherenkovFreeDirection >
    >
    {
        using LowerMargin = typename pmacc::math::CT::make_Int<
            simDim,
            1
        >::type;
        using UpperMargin = typename pmacc::math::CT::make_Int<
            simDim,
            2
        >::type;
    };

} // namespace traits
} // namespace picongpu
