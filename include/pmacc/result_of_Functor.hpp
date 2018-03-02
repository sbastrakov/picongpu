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


namespace pmacc
{
namespace result_of
{

template<
    typename _Functor,
    typename Arg0 = boost::mp11::mp_void<>,
    typename Arg1 = boost::mp11::mp_void<>,
    typename Arg2 = boost::mp11::mp_void<>,
    typename Arg3 = boost::mp11::mp_void<>,
    typename Arg4 = boost::mp11::mp_void<>,
    typename Arg5 = boost::mp11::mp_void<>,
    typename Arg6 = boost::mp11::mp_void<>,
    typename Arg7 = boost::mp11::mp_void<>,
    typename Arg8 = boost::mp11::mp_void<>,
    typename Arg9 = boost::mp11::mp_void<>,
    typename Arg10 = boost::mp11::mp_void<>,
    typename Arg11 = boost::mp11::mp_void<>,
    typename Arg12 = boost::mp11::mp_void<>,
    typename dummy = boost::mp11::mp_void<>
>
struct Functor
{
    typedef typename _Functor::result_type type;
};

} // result_of
} // PMacc
