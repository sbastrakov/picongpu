/* Copyright 2014-2018 Rene Widera, Benjamin Worpitz
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

#include "pmacc/compileTime/conversion/SeqToMap.hpp"
#include "pmacc/compileTime/conversion/TypeToAliasPair.hpp"
#include "pmacc/compileTime/conversion/TypeToPair.hpp"
#include "pmacc/compileTime/conversion/MakeSeqFromNestedSeq.hpp"
#include "pmacc/compileTime/conversion/MakeSeq.hpp"
#include "pmacc/math/Vector.hpp"
#include "pmacc/types.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/bind.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>
#include <boost/type_traits/is_same.hpp>

#include <type_traits>


namespace pmacc
{

namespace detail
{
/** Create tuples out of the elements of N sequences
 *
 * Combines all elements of N given sequences in T_MplSeq into N-tuples.
 * If the number of elements in each sequence is S0, S1, ... S(N-1)
 * than the resulting sequence will contain S0 * S1 * ... S(N-1) tuples.
 *
 * @tparam T_MplSeq sequence of input sequences
 * @tparam T_TmpResult temporary result
 * @tparam T_isEmpty true if T_MplSeq is empty else false
 */
template<typename T_MplSeq,
typename T_TmpResult = bmp11::mp_list<>,
bool T_isEmpty = bmp11::mp_empty<T_MplSeq>::value
>
struct AllCombinations;

/** implementation for inner recursive creation
 */
template<typename T_MplSeq, typename T_TmpResult>
struct AllCombinations<T_MplSeq, T_TmpResult, false >
{
    typedef T_MplSeq MplSeq;
    typedef T_TmpResult TmpResult;

    static constexpr uint32_t rangeVectorSize = bmp11::mp_size<MplSeq>::value;
    using LastElement = bmp11::mp_at_c<MplSeq, rangeVectorSize - 1 >;
    typedef bmp11::mp_empty<LastElement> IsLastElementEmpty;
    typedef typename MakeSeq<LastElement>::type LastElementAsSequence;
    using ShrinkedRangeVector = bmp11::mp_take_c<
        MplSeq,
        rangeVectorSize - 1
    >;

    /** Assign to each element in a sequence of CT::Vector(s) a type at a given
    *  component position
    *
    * @tparam T_ComponentPos position of the component to be changed (type must be std::integral_constant<uint32_t,X>)
    * @tparam T_Element value (type) which should replace the component at position T_Component
    *                   in the CT::Vector elements
    */
    template<
        typename T_ComponentPos,
        typename T_Element
    >
    struct AssignToAnyElementInVector
    {
        typedef TmpResult InVector;
        typedef T_Element Element;

        using type = bmp11::mp_transform<
            pmacc::math::CT::Assign<
                bmp11::_1,
                T_ComponentPos,
                Element
            >,
            InVector
        >;
    };

    using NestedSeq = bmp11::mp_transform<
        AssignToAnyElementInVector<
            std::integral_constant<uint32_t, rangeVectorSize - 1 >,
            bmp11::_1
        >,
        LastElementAsSequence
    >;

    typedef typename MakeSeqFromNestedSeq<NestedSeq>::type OneSeq;

    typedef typename detail::AllCombinations<ShrinkedRangeVector, OneSeq>::type ResultIfNotEmpty;
    typedef typename bmp11::mp_if<IsLastElementEmpty,bmp11::mp_list<>,ResultIfNotEmpty> type;
};

/** recursive end implementation
 */
template<typename T_MplSeq, typename T_TmpResult>
struct AllCombinations<T_MplSeq, T_TmpResult, true >
{
    typedef T_TmpResult type;
};

} //detail


/** Create tuples out of the elements of N sequences
 *
 * Combines all elements of N given sequences in T_MplSeq into N-tuples.
 * If the number of elements in each sequence is S0, S1, ... S(N-1)
 * than the resulting sequence will contain S0 * S1 * ... S(N-1) tuples.
 *
 * example:
 *
 * sequence  == [ ]
 * tuple     == ( )
 *
 * T_MplSeq = [[1,2],[1],[4,3]]
 * combined to
 * AllCombinations<T_MplSeq>::type = [(1,1,4),(1,1,3),(2,1,4),(2,1,3)]
 *
 * @tparam T_MplSeq N-dimensional sequence with input values
 *                  or single type (e.g. `std::integral_constant<uint32_t,5>`)
 *                  (if `T_MplSeq` is only one type it will be transformed to a sequence)
 * @typedef AllCombinations<T_MplSeq>::type
 *          MplSequence of N-tuples
 */
template<typename T_MplSeq>
struct AllCombinations
{
    /* if T_MplSeq is no sequence it is a single type, we put this type in
     * a sequence because all next algorithms can only work with sequences */
    typedef typename MakeSeq<T_MplSeq>::type MplSeq;

    static constexpr uint32_t rangeVectorSize = bmp11::mp_size<MplSeq>::value;
    using LastElement = bmp11::mp_at_c<MplSeq, rangeVectorSize - 1 >;
    typedef bmp11::mp_empty<LastElement> IsLastElementEmpty;
    typedef typename MakeSeq<LastElement>::type LastElementAsSequence;

    using ShrinkedRangeVector = bmp11::mp_take_c<
        MplSeq
        rangeVectorSize - 1
    >;

    /* transform all elements in the vector to math::CT::vector<> */
    typedef math::CT::Vector<> EmptyVector;
    using FirstList = bmp11::mp_transform<
        pmacc::math::CT::Assign<EmptyVector, std::integral_constant<uint32_t, rangeVectorSize - 1 >, bmp11::_1>,
        LastElementAsSequence
    >;

    /* result type: MplSequence of N-tuples */
    typedef typename detail::AllCombinations<ShrinkedRangeVector, FirstList>::type ResultIfNotEmpty;
    typedef typename bmp11::mp_if<IsLastElementEmpty,bmp11::mp_list<>,ResultIfNotEmpty> type;
};


}//namespace pmacc
