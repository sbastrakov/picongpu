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
 * \tparam T_List list of keys to search in
 * \tparam T_Key key or alias of a key in the sequence
 * \tparam T_KeyNotFoundPolicy metafunction to be called if key is not found
 */
template<
    typename T_List,
    typename T_Key,
    template<
        typename T_List,
        typename T_Key
    > class T_KeyNotFoundPolicy = errorHandlerPolicies::ReturnType
>
struct GetKeyFromAlias
{
private:

    /*create a map where Key is a undeclared alias and value is real type*/
    using AliasMap = SeqToMap<
        T_List,
        TypeToAliasPair_t
    >;

    /*create a map where Key and value is real type*/
    using KeyMap = SeqToMap<
        T_List,
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
        T_KeyNotFoundPolicy<T_List, T_Key>,
        MapType
    >;
};

template<
    typename T_List,
    typename T_Key,
    template<
        typename T_List,
        typename T_Key
    > class T_KeyNotFoundPolicy = errorHandlerPolicies::ReturnType
>
using GetKeyFromAlias_t = typename GetKeyFromAlias<
    T_List,
    T_Key,
    T_KeyNotFoundPolicy
>::type;

}//namespace pmacc
