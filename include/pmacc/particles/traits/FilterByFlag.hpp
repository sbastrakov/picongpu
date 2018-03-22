/* Copyright 2015-2018 Heiko Burau
 *
 * This file is part of PMacc.
 *
 * PMacc is free software: you can redistribute it and/or modify
 * it under the terms of either the GNU General Public License or
 * the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PMacc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License and the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and the GNU Lesser General Public License along with PMacc.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once


#include "pmacc/types.hpp"
#include "pmacc/traits/HasFlag.hpp"

#include <boost/mp11/algorithm.hpp>


namespace pmacc
{
namespace particles
{
namespace traits
{
namespace detail
{

template<
    typename T_List,
    typename T_Flag
>
struct FilterByFlag
{
    template< typename T_Species >
    using HasFlag = typename ::pmacc::traits::HasFlag<
        typename T_Species::FrameType,
        T_Flag
    >::type;

    using type = bmp11::mp_copy_if<
        T_List,
        HasFlag
    >;
};

} // namespace detail

/** Return a new list of particle species carrying the flag.
 *
 * @tparam T_List list of particle species
 * @tparam T_Flag flag to be filtered
 */
template<
    typename T_List,
    typename T_Flag
>
using FilterByFlag = typename detail::FilterByFlag<
    T_List,
    T_Flag
>::type;

} //namespace traits
} //namespace particles
} //namespace pmacc
