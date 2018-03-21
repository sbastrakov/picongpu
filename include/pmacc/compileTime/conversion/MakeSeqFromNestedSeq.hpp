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


#include "pmacc/compileTime/conversion/ToSeq.hpp"
#include "pmacc/compileTime/conversion/JoinToSeq.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>


namespace pmacc
{

/** Combine all elements of the input type or list to a single list
 *
 * If the input is a list and its elements are lists themselves, all of their
 * elements will be added to the resulting list
 */
template< typename T >
using MakeSeqFromNestedSeq = bmp11::mp_fold<
    ToSeq< T >,
    bmp11::mp_list<>,
    JoinToSeq
>;

} //namespace pmacc
