/* Copyright 2013-2020 Axel Huebl, Heiko Burau, Rene Widera
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

#include "picongpu/algorithms/DifferenceToUpper.hpp"
#include "picongpu/algorithms/DifferenceToLower.hpp"
#include "picongpu/fields/MaxwellSolver/Yee/Curl.def"

#include <pmacc/types.hpp>


namespace picongpu
{
namespace fields
{
namespace maxwellSolver
{
namespace yee
{

    template< typename T_Difference >
    struct Curl
    {
        using Difference = T_Difference;
        using LowerMargin = typename Difference::OffsetOrigin;
        using UpperMargin = typename Difference::OffsetEnd;

        struct Components
        {
            float_X dFxDy;
            float_X dFxDz;
            float_X dFyDx;
            float_X dFyDz;
            float_X dFzDx;
            float_X dFzDy;
        };

        template<class Memory >
        HDINLINE Components getComponents( const Memory & mem ) const
        {
            using DifferentiatorX = typename Difference::template GetDifference< 0 >;
            using DifferentiatorY = typename Difference::template GetDifference< 1 >;
            using DifferentiatorZ = typename Difference::template GetDifference< 2 >;
            auto const dFDx = DifferentiatorX{}( mem );
            auto const dFDy = DifferentiatorY{}( mem );
            auto const dFDz = DifferentiatorZ{}( mem );
            Components result;
            result.dFzDy = dFDy.z();
            result.dFyDz = dFDz.y();
            result.dFxDz = dFDz.x();
            result.dFzDx = dFDx.z();
            result.dFyDx = dFDx.y();
            result.dFxDy = dFDy.x();
            return result;
        }

        template<class Memory >
        HDINLINE typename Memory::ValueType operator()( Memory const & mem ) const
        {
            auto const components = getComponents( mem );
            return float3_X(
                components.dFyDz - components.dFzDy,
                components.dFzDx - components.dFxDz,
                components.dFxDy - components.dFyDx
            );
        }
    };

} // namespace yee
} // namespace maxwellSolver
} // namespace fields
} // namespace picongpu
