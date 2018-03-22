/* Copyright 2013-2018 Rene Widera, Felix Schmitt, Alexander Grund
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

#include "pmacc/compileTime/GetKeyFromAlias.hpp"
#include "pmacc/compileTime/errorHandlerPolicies/ThrowValueNotFound.hpp"

#include <boost/mp11/algorithm.hpp>


namespace pmacc
{

namespace detail
{

template<
    typename T_List,
    typename T_Lookup,
    template<
        typename T_List,
        typename T_Key
    > class T_AliasNotFoundPolicy = errorHandlerPolicies::ThrowValueNotFound
>
struct ResolveAliases
{
    template< typename T_Identifier >
    using GetKeyFromAliasAccessor = GetKeyFromAlias_t<
        T_Lookup,
        T_Identifier,
        T_AliasNotFoundPolicy
    >;

    using type = bmp11::mp_transform<
        GetKeyFromAliasAccessor,
        T_List
    >;
};

} // namespace detail

/** Translate all pmacc alias types to full specialized types
 *
 * Use lookup list to translate types
 * The policy is used if the type from T_List is not in T_Lookup a compile time error is triggered
 *
 * @tparam T_List source list with types to translate
 * @tparam T_Lookup lookup list to translate aliases
 */
template<
    typename T_List,
    typename T_Lookup,
    template<
        typename T_List,
        typename T_Key
    > class T_AliasNotFoundPolicy = errorHandlerPolicies::ThrowValueNotFound
>
using ResolveAliases = typename detail::ResolveAliases<
    T_List,
    T_Lookup,
    T_AliasNotFoundPolicy
>::type;

} //namespace pmacc
