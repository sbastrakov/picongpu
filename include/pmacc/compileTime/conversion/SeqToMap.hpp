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
#include "pmacc/compileTime/accessors/Identity.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/map.hpp>


namespace pmacc
{

namespace detail
{

template<
    typename T_List,
    template< typename > class T_UnaryOperator,
    template< typename > class T_Accessor
>
struct SeqToMap
{
    template< typename T >
    using Op = T_UnaryOperator< T_Accessor< T > >;

    using Lists = bmp11::mp_transform<
        Op,
        T_List
    >;

    using type = bmp11::mp_fold<
        Lists,
        bmp11::mp_list<>,
        bmp11::mp_map_insert
    >;
};

} // namespace detail

/** convert boost mp11 list to an mp11 map
*
* @tparam T_List boost mp11 list
* @tparam T_UnaryOperator unary metafunction to translate a type from T_List to an mp11 list
* @tparam T_Accessor an unary metafunction which is used before the type
* from the list is passed to T_UnaryOperator
*/
template<
    typename T_List,
    template< typename > class T_UnaryOperator,
    template< typename > class T_Accessor = compileTime::accessors::Identity_t
>
using SeqToMap = typename detail::SeqToMap<
    T_List,
    T_UnaryOperator,
    T_Accessor
>::type;

} //namespace pmacc
