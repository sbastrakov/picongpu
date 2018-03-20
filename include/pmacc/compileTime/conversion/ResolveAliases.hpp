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
#include <boost/mp11/utility.hpp>


namespace pmacc
{

/** Translate all pmacc alias types to full specialized types
 *
 * Use lookup list to translate types
 * The policy is used if the type from T_Seq is not in T_SeqLookup a compile time error is triggered
 *
 * @tparam T_Seq source list with types to translate
 * @tparam T_SeqLookup lookup list to translate aliases
 */
template<
    typename T_Seq,
    typename T_SeqLookup,
    typename T_AliasNotFoundPolicy = errorHandlerPolicies::ThrowValueNotFound
>
struct ResolveAliases
{
    typedef T_Seq Seq;
    typedef T_SeqLookup SeqLookup;
    typedef T_AliasNotFoundPolicy AliasNotFoundPolicy;

    template< typename T_Identifier >
    using GetKeyFromAliasAccessor = GetKeyFromAlias_t< SeqLookup, T_Identifier, AliasNotFoundPolicy >;

    using type = bmp11::mp_transform<
        GetKeyFromAliasAccessor,
        Seq
    >;
};

}//namespace pmacc
