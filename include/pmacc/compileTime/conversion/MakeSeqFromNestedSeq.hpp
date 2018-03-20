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
#include <boost/mp11/utility.hpp>


namespace pmacc
{

/** combine all elements of the input type to a single vector
 *
 * If elements of the input sequence are a sequence themself, all of their
 * elements will be added to the resulting sequence
 *
 * @tparam T_In a boost mpl sequence or single type
 */
template<typename T_In>
struct MakeSeqFromNestedSeq
{
public:
    using type = bmp11::mp_fold<
      ToSeq<T_In >,
      bmp11::mp_list<>,
      JoinToSeq
    >;
};

} //namespace pmacc
