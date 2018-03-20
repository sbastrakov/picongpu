/* Copyright 2013-2018 Rene Widera, Benjamin Worpitz, Alexander Grund
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
#include "pmacc/compileTime/conversion/SeqToMap.hpp"
#include "pmacc/compileTime/conversion/TypeToAliasPair.hpp"
#include "pmacc/compileTime/conversion/TypeToPair.hpp"
#include "pmacc/compileTime/errorHandlerPolicies/ReturnType.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/map.hpp>
#include <boost/mp11/utility.hpp>


namespace pmacc
{

/**
 * Returns the key type from an alias
 *
 * \tparam T_MPLSeq Sequence of keys to search
 * \tparam T_Key Key or alias of a key in the sequence
 * \tparam T_KeyNotFoundPolicy Binary meta-function that is called like (T_MPLSeq, T_Key)
 *         when T_Key is not found in the sequence. Default is to return bmp11::mp_void<>
 */
template<typename T_MPLSeq,
         typename T_Key,
         typename T_KeyNotFoundPolicy = errorHandlerPolicies::ReturnType<void, T_MPLSeq, T_Key>
>
struct GetKeyFromAlias
{
private:

    /*create a map where Key is a undeclared alias and value is real type*/
    using AliasMap = SeqToMap_t<
        T_MPLSeq,
        TypeToAliasPair_t
    >;

    /*create a map where Key and value is real type*/
    using KeyMap = SeqToMap_t<
        T_MPLSeq,
        TypeToPair
    >;

    /*combine both maps*/
    using FullMap = bmp11::mp_fold<
        AliasMap,
        KeyMap,
        bmp11::mp_map_insert
    >;

    /* search for given key,
     * - we get the real type if key found
     * - else we get boost::mp11::mp_void<>
     */
    using MapType = bmp11::mp_map_find<
        FullMap,
        T_Key
    >;

public:

    using type = bmp11::mp_if<
        std::is_same<
            MapType,
            void
        >,
        T_KeyNotFoundPolicy,
        MapType
    >;
};

template<
    typename T_Seq,
    typename T_Key,
    typename T_KeyNotFoundPolicy = errorHandlerPolicies::ReturnType<>
>
using GetKeyFromAlias_t = typename GetKeyFromAlias<
    T_Seq,
    T_Key,
    T_KeyNotFoundPolicy
>::type;

}//namespace pmacc
