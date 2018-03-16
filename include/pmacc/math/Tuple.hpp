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

#include "pmacc/types.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_shifted_params.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/static_assert.hpp>

#include <type_traits>


namespace pmacc
{
namespace math
{

#ifndef TUPLE_MAX_DIM
#define TUPLE_MAX_DIM 8
#endif

#define CONSTRUCTOR(Z, N, _)                                \
    template<BOOST_PP_ENUM_PARAMS(N, typename Arg)>         \
    HDINLINE                                                \
    Tuple(BOOST_PP_ENUM_BINARY_PARAMS(N, const Arg, &arg))  \
    : value(arg0),                                          \
      base(BOOST_PP_ENUM_SHIFTED_PARAMS(N, arg))            \
    {                                                       \
        BOOST_STATIC_ASSERT(dim == N);                      \
    }


template<typename TypeList, bool ListEmpty = bmp11::mp_empty<TypeList>::value>
class Tuple;

template<typename TypeList>
class Tuple<TypeList, true> {};

template<typename TypeList>
class Tuple<TypeList, false>
    : public Tuple<bmp11::mp_pop_front<TypeList>>
{
public:
    static constexpr int dim = bmp11::mp_size<TypeList>::value;
    typedef TypeList TypeList_;
private:
    typedef Tuple<bmp11::mp_pop_front<TypeList>> base;

    typedef typename bmp11::mp_front<TypeList> Value;
    typedef typename boost::remove_reference<Value>::type pureValue;

    Value value;
public:
    HDINLINE Tuple() {}

    HDINLINE Tuple(Value arg0) : value(arg0)
    {
        BOOST_STATIC_ASSERT(dim == 1);
    }

    BOOST_PP_REPEAT_FROM_TO(2, BOOST_PP_INC(TUPLE_MAX_DIM), CONSTRUCTOR, _)

    template<int i>
    HDINLINE
    bmp11::mp_at_c<TypeList, i>&
    at_c()
    {
        return this->at(bmp11::mp_int<i>::value);
    }
    template<int i>
    HDINLINE
    const bmp11::mp_at_c<TypeList, i>&
    at_c() const
    {
        return this->at(bmp11::mp_int<i>::value);
    }

    HDINLINE Value& at(bmp11::mp_int<0>)
    {
        return value;
    }

    HDINLINE const Value& at(bmp11::mp_int<0>) const
    {
        return value;
    }

    template<typename Idx>
    HDINLINE
    bmp11::mp_at<TypeList, Idx>&
    at(Idx)
    {
        return base::at(bmp11::mp_int<Idx::value - 1>);
    }

    template<typename Idx>
    HDINLINE
    const bmp11::mp_at<TypeList, Idx>&
    at(Idx) const
    {
        return base::at(bmp11::mp_int<Idx::value - 1>);
    }
};

#undef CONSTRUCTOR

#define MAKE_TUPLE(Z, N, _) \
    template<BOOST_PP_ENUM_PARAMS(N, typename Value)> \
    HDINLINE \
    Tuple<bmp11::mp_list<BOOST_PP_ENUM_PARAMS(N, Value)> > \
    make_Tuple(BOOST_PP_ENUM_BINARY_PARAMS(N, Value, value)) \
    { \
        return Tuple<bmp11::mp_list<BOOST_PP_ENUM_PARAMS(N, Value)> > \
            (BOOST_PP_ENUM_PARAMS(N, value)); \
    }

BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_INC(TUPLE_MAX_DIM), MAKE_TUPLE, _)

#undef MAKE_TUPLE

namespace result_of
{
template<typename TTuple, int i>
struct at_c
{
    using type = bmp11::mp_at_c<typename TTuple::TypeList_, i>;
};
} // result_of

} // math
} // PMacc
