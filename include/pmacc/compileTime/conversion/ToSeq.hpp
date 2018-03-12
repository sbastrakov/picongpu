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

#include "pmacc/types.hpp"

#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>


namespace pmacc
{

/** cast type to boost mp11 mp_list
 * @return ::type if T_Type is sequence then identity of T_Type
 *                else boost::mp11::mp_list<T_Type>
 */
template<typename T_Type>
struct ToSeq
{
    using type = bmp11::mp_if<
        bmp11::mp_is_list< T_Type >,
        T_Type,
        bmp11::mp_list<T_Type>
    >;
};

}//namespace pmacc
