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

#include <boost/mpl/front.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>


namespace pmacc
{
namespace particles
{
namespace compileTime
{

    /* find a type within a sequence by name or the type itself
     *
     * pmacc::traits::GetCTName is used to translate each element of
     * T_MPLSeq into a name.
     *
     * @tparam T_MPLSeq source sequence where we search T_Identifier
     * @tparam T_Identifier name or type to search
     */
    template<
        typename T_MPLSeq,
        typename T_Identifier,
        typename T_KeyNotFoundPolicy = pmacc::errorHandlerPolicies::ThrowValueNotFound
    >
    struct FindByNameOrType
    {
        using KeyNotFoundPolicy = T_KeyNotFoundPolicy;

        template< typename T_Value >
        struct HasTypeOrName
        {
            using type = bmp11::mp_or<
                bmp11::mp_same<
                    T_Identifier,
                    T_Value
                >,
                bmp11::mp_same<
                    pmacc::traits::GetCTName_t< T_Value >,
                    T_Identifier
                >
            >;
        };

        using FilteredSeq = bmp11::mp_copy_if<
            T_MPLSeq,
            bmp11::mp_identity_t< HasTypeOrName >
        >;

        using type = typename bmp11::mp_if<
            bmp11::mp_empty< FilteredSeq >,
            KeyNotFoundPolicy<
                T_MPLSeq,
                T_Identifier
            >,
            bmpl::front< FilteredSeq >
        >::type;
    };

    template<
        typename T_MPLSeq,
        typename T_Identifier,
        typename T_KeyNotFoundPolicy = pmacc::errorHandlerPolicies::ThrowValueNotFound
    >
    using FindByNameOrType_t = typename FindByNameOrType<
        T_MPLSeq,
        T_Identifier,
        T_KeyNotFoundPolicy
    >::type;

} // namespace compileTime
} // namespace particles
} // namespace pmacc
