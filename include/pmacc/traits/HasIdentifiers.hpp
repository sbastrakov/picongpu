/* Copyright 2017-2018 Axel Huebl
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

#include "pmacc/traits/HasIdentifier.hpp"
#include <boost/mp11/bind.hpp>


namespace pmacc
{
namespace traits
{

    /** Checks if an object has all specified identifiers
     *
     * @tparam T_Object any object (class or typename)
     * @tparam T_SeqKeys a sequence of identifiers
     *
     * This struct must define
     * ::type (boost::mp11::mp_bool<>)
     */
    template<
        typename T_Object,
        typename T_SeqKeys
    >
    struct HasIdentifiers
    {
        using type = bmp11::mp_all_of_q<
            T_SeqKeys,
            bmp11::mp_bind_front<
                HasIdentifier_t,
                T_Object
            >
        >;
    };

    template<
        typename T_Object,
        typename T_SeqKeys
    >
    bool hasIdentifiers(
        T_Object const &,
        T_SeqKeys const &
    )
    {
        return HasIdentifiers<
            T_Object,
            T_SeqKeys
        >::type::value;
    }

} // namespace traits
} // namespace pmacc
