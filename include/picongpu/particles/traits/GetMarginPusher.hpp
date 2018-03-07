/* Copyright 2015-2018 Richard Pausch
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
#include "picongpu/traits/GetMargin.hpp"
#include "picongpu/particles/traits/GetInterpolation.hpp"
#include "picongpu/particles/traits/GetPusher.hpp"

#include <boost/mp11/bind.hpp>


namespace picongpu
{

namespace traits
{
template<typename T_Species>
struct GetMarginPusher
{
    typedef pmacc::math::CT::add<
        GetLowerMargin< GetInterpolation< bmp11::_1 > >,
        GetLowerMargin< GetPusher< bmp11::_1 > >
    > AddLowerMargins;
    typedef typename bmpl::apply<AddLowerMargins, T_Species>::type LowerMargin;

    typedef pmacc::math::CT::add<
        GetUpperMargin< GetInterpolation< bmp11::_1 > >,
        GetUpperMargin< GetPusher< bmp11::_1 > >
    > AddUpperMargins;
    typedef typename bmpl::apply<AddUpperMargins, T_Species>::type UpperMargin;
};

template<typename T_Species>
struct GetLowerMarginPusher
{
    using type = typename traits::GetMarginPusher<T_Species>::LowerMargin;
};

template<typename T_Species>
struct GetUpperMarginPusher
{
    using type = typename traits::GetMarginPusher<T_Species>::UpperMargin;
};

}// namespace traits
}// namespace picongpu
