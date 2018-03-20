/* Copyright 2013-2018 Rene Widera
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

#include <boost/mp11/list.hpp>

#include "pmacc/compileTime/conversion/ToSeq.hpp"


namespace pmacc
{

/** Join both input types to one boost mp11 list
 *
 * @tparam T_1 a boost mp11 list or single type
 * @tparam T_2 a boost mp11 list or single type
 */
template<
    typename T_1,
    typename T_2 = bmp11::mp_list<>
>
using JoinToSeq = bmp11::mp_append<
    ToSeq< T_1 >,
    ToSeq< T_2 >
>;

} //namespace pmacc
