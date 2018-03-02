/* Copyright 2013-2018 Heiko Burau, Rene Widera
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

#include <boost/mp11/function.hpp>

#include <stdint.h>


namespace pmacc
{
namespace math
{
namespace CT
{

template<
    typename X = mp11::mp_void,
    typename Y = mp11::mp_void,
    typename Z = mp11::mp_void
>
struct Float
{
    typedef X x;
    typedef Y y;
    typedef Z z;

    static constexpr int dim = 3;
};

template<>
struct Float<> {};

template< typename X >
struct Float< X >
{
    typedef X x;

    static constexpr int dim = 1;
};

template<
    typename X,
    typename Y
>
struct Float<
    X,
    Y
>
{
    typedef X x;
    typedef Y y;

    static constexpr int dim = 2;
};

} // CT
} // math
} // pmacc
