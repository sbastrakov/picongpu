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


#include "pmacc/compileTime/accessors/Identity.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/utility.hpp>


namespace pmacc
{
namespace detail
{

    struct Empty {};

    /** get combined type which inherit from a boost mp11 list
     *
     * @tparam T_Sequence boost mp11 list with classes
     * @tparam T_Accessor unary operator to transform each element of the sequence
     */
    template<
        typename T_Sequence,
        template< typename > class T_Accessor = compileTime::accessors::Identity
    >
    using InheritLinearly = bmp11::mp_fold<
        bmp11::mp_transform<
            T_Sequence,
            T_Accessor
        >
        Empty,
        bmp11::mp_inherit
    >;

} //namespace detail

    /** type which inherits from multiple classes
     *
     * @tparam T_Sequence boost mpl sequence with classes
     * @tparam T_Accessor unary operator to transform each element of the sequence
     */
    template<
        typename T_Sequence,
        template< typename > class T_Accessor = compileTime::accessors::Identity
    >
    struct InheritLinearly : detail::InheritLinearly<
        T_Sequence,
        T_Accessor
    >
    {
    };

} //namespace pmacc
