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


#include "pmacc/compileTime/conversion/MakeSeqFromNestedSeq.hpp"

#include <boost/mp11/list.hpp>


namespace pmacc
{

/** combine all input types to one sequence
 *
 * Note: if the input type is a sequence itself, its elements will be unfolded
 *       and added separately
 *
 * @tparam T_Args a boost mpl sequence or single type
 *
 * @code
 * using MyType = typename MakeSeq< A, B >::type
 * using MyType2 = typename MakeSeq< boost::mp11::mp_list<A, B>, C >::type
 * @endcode
 *
 */
template< typename... T_Args >
struct MakeSeq
{
    typedef typename MakeSeqFromNestedSeq<
        bmp11::mp_list< T_Args... >
    >::type type;
};

/** short hand definition for @see MakeSeq<> */
template< typename... T_Args >
using MakeSeq_t = typename MakeSeq< T_Args... >::type;

} //namespace pmacc
