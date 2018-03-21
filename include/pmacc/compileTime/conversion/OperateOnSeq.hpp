/* Copyright 2015-2018 Rene Widera
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

#include "pmacc/compileTime/accessors/Identity.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/utility.hpp>


namespace pmacc
{

/** Transform boost mp11 list into another list by applying operations
 *
 * @tparam T_List boost mp11 list
 * @tparam T_UnaryOperator unary operator to apply for each element
 * @tparam T_Accessor an unary lambda operator that is used before the type
 * from the sequence is passed to T_UnaryOperator
 */
template<
    typename T_List,
    typename T_UnaryOperator,
    typename T_Accessor = compileTime::accessors::Identity_t
>
using OperateOnSeq = bmp11::mp_transform<
    T_UnaryOperator,
    bmp11::mp_transform<
        T_Accessor,
        T_List
    >
>;

}//namespace pmacc
