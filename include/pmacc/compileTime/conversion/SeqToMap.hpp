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
#include <boost/mp11/bind.hpp>
#include <boost/mp11/map.hpp>
#include <boost/mp11/utility.hpp>


namespace pmacc
{

/** convert boost mp11 list to an mp11 map
 *
 * @tparam T_MPLSeq boost mp11 list
 * @tparam T_UnaryOperator unary operator to translate type from the sequence
 * to a mp11 list
 * @tparam T_Accessor An unary lambda operator which is used before the type
 * from the sequence is passed to T_UnaryOperator
 * @return ::type mp11 map
 */
template<
    typename T_MPLSeq,
    template< typename > class T_UnaryOperator,
    typename T_Accessor = compileTime::accessors::Identity<>
>
struct SeqToMap
{
    template< typename T >
    using Op = T_UnaryOperator< typename T_Accessor< X >::type >;

    using Lists = bmp11::mp_transform<
        Op,
        T_MPLSeq
    >;
    using type = bmp11::mp_fold<
        Lists,
        bmp11::mp_list<>,
        bmp11::mp_map_insert
    >;
};

template<
    typename T_MPLSeq,
    template< typename > class T_UnaryOperator,
    typename T_Accessor = compileTime::accessors::Identity<>
>
using SeqToMap_t = typename SeqToMap< T_MPLSeq, T_UnaryOperator, T_Accessor >::type;

}//namespace pmacc
