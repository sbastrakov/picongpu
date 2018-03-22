/* Copyright 2018 Rene Widera
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

#include "pmacc/traits/GetCTName.hpp"
#include "pmacc/compileTime/errorHandlerPolicies/ThrowValueNotFound.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/function.hpp>
#include <boost/mp11/list.hpp>


namespace pmacc
{
namespace particles
{
namespace compileTime
{
namespace detail
{

template<
    typename T_List,
    typename T_Identifier,
    template<
        typename T_List,
        typename T_Identifier
    > class T_KeyNotFoundPolicy
>
struct FindByNameOrType
{
    template< typename T_Value >
    using HasTypeOrName = bmp11::mp_or<
        bmp11::mp_same<
            T_Identifier,
            T_Value
        >,
        bmp11::mp_same<
            pmacc::traits::GetCTName_t< T_Value >,
            T_Identifier
        >
    >;
 
    using FilteredSeq = bmp11::mp_copy_if<
        T_List,
        HasTypeOrName
    >;

    using type = bmp11::mp_if<
        bmp11::mp_empty< FilteredSeq >,
        T_KeyNotFoundPolicy<
            T_List,
            T_Identifier
        >,
        bmp11::mp_front< FilteredSeq >
    >;
};

} // namespace detail

/**
 * Find a type within a sequence by name or the type itself
 *
 * pmacc::traits::GetCTName is used to translate each element of
 * T_List into a name.
 *
 * @tparam T_List source list where we search T_Identifier
 * @tparam T_Identifier name or type to search
 */
template<
    typename T_List,
    typename T_Identifier,
    template<
        typename T_List,
        typename T_Identifier
    > class T_KeyNotFoundPolicy = pmacc::errorHandlerPolicies::ThrowValueNotFound
>
using FindByNameOrType = typename detail::FindByNameOrType<
    T_List,
    T_Identifier,
    T_KeyNotFoundPolicy
>::type;

} // namespace compileTime
} // namespace particles
} // namespace pmacc
