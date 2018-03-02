/* Copyright 2013-2018 Heiko Burau, Rene Widera, Benjamin Worpitz
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

#include "pmacc/math/Vector.hpp"
#include "pmacc/types.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/function.hpp>
#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>

#include <stdint.h>
#include <type_traits>

namespace pmacc
{
namespace math
{
namespace CT
{

namespace mp11 = boost::mp11;

namespace detail
{
template<int dim>
struct VectorFromCT;

template<>
struct VectorFromCT<1>
{

    template<typename Vec, typename CTVec>
    HDINLINE void operator()(Vec& vec, CTVec) const
    {
        BOOST_STATIC_ASSERT(Vec::dim == 1);
        BOOST_STATIC_ASSERT(CTVec::dim == 1);
        vec[0] = (typename Vec::type)CTVec::x::value;
    }
};

template<>
struct VectorFromCT<2>
{

    template<typename Vec, typename CTVec>
    HDINLINE void operator()(Vec& vec, CTVec) const
    {
        BOOST_STATIC_ASSERT(Vec::dim == 2);
        BOOST_STATIC_ASSERT(CTVec::dim == 2);
        vec[0] = (typename Vec::type)CTVec::x::value;
        vec[1] = (typename Vec::type)CTVec::y::value;
    }
};

template<>
struct VectorFromCT<3>
{

    template<typename Vec, typename CTVec>
    HDINLINE void operator()(Vec& vec, CTVec) const
    {
        BOOST_STATIC_ASSERT(Vec::dim == 3);
        BOOST_STATIC_ASSERT(CTVec::dim == 3);
        vec[0] = (typename Vec::type)CTVec::x::value;
        vec[1] = (typename Vec::type)CTVec::y::value;
        vec[2] = (typename Vec::type)CTVec::z::value;
    }
};

struct VectorParameterPlaceholder {};

template<typename Arg0>
struct TypeSelector
{
    typedef Arg0 type;
};

/** get integral type*/
template<typename T, T value>
struct TypeSelector<std::integral_constant<T, value > >
{
    typedef T type;
};

template<>
struct TypeSelector<VectorParameterPlaceholder>
{
    typedef mp11::mp_int<0> type;
};

}

namespace mp11 = boost::mp11;


template<
    typename Arg0 = detail::VectorParameterPlaceholder,
    typename Arg1 = detail::VectorParameterPlaceholder,
    typename Arg2 = detail::VectorParameterPlaceholder
>
struct Vector
{
    typedef Arg0 x;
    typedef Arg1 y;
    typedef Arg2 z;

    using mplVector = mp11::mp_if<
        mp11::mp_same<
            x,
            detail::VectorParameterPlaceholder
        >,
        mp11::mp_list<>,
        mp11::mp_if<
            mp11::mp_same<
                y,
                detail::VectorParameterPlaceholder
            >,
            mp11::mp_list< x >,
            mp11::mp_if<
                mp11::mp_same <
                    z,
                    detail::VectorParameterPlaceholder
                >,
                mp11::mp_list<
                    x,
                    y
                >,
                mp11::mp_list<
                    x,
                    y,
                    z
                >
            >
        >
    >;

    template<int element>
    struct at
    {
        using type = mp11::mp_at_c<mplVector, element>;
    };

    static constexpr int dim = mp11::mp_size<mplVector >::value;

    typedef typename detail::TypeSelector<x>::type type;
    typedef Vector<x, y, z> This;
    typedef math::Vector<type, dim> RT_type;
    typedef This vector_type;

    template<typename OtherType>
    HDINLINE
    operator math::Vector<OtherType, dim>() const
    {
        math::Vector<OtherType, dim> result;
        math::CT::detail::VectorFromCT<dim>()(result, *this);
        return result;
    }

    /** Create a runtime Vector
     *
     *  Creates the corresponding runtime vector object.
     *
     *  \return RT_type runtime vector with same value type
     */
    static HDINLINE RT_type toRT()
    {
        math::Vector<type, dim> result;
        math::CT::detail::VectorFromCT<dim>()(result, This());
        return result;
    }
};

//*********************************************************

//________________________OperatorBase____________________________

template<
    typename Lhs,
    typename Rhs,
    typename T_BinaryOperator
>
struct applyOperator
{
    using type = typename applyOperator<
        typename Lhs::vector_type,
        typename Rhs::vector_type,
        T_BinaryOperator
    >::type;
};

template<
    typename T_TypeA,
    typename T_TypeB,
    typename T_BinaryOperator
>
struct applyOperator<
    CT::Vector< T_TypeA >,
    CT::Vector< T_TypeB >,
    T_BinaryOperator
>
{
    using OpResult = T_BinaryOperator<
        T_TypeA,
        T_TypeB
    >;
    using type = CT::Vector< OpResult >;
};

template<
    typename T_TypeA0,
    typename T_TypeA1,
    typename T_TypeB0,
    typename T_TypeB1,
    typename T_BinaryOperator
>
struct applyOperator<
    CT::Vector<
        T_TypeA0,
        T_TypeA1
    >,
    CT::Vector<
        T_TypeB0,
        T_TypeB1
    >,
    T_BinaryOperator
>
{
    using OpResult0 = T_BinaryOperator<
        T_TypeA0,
        T_TypeB0
    >;
    using OpResult1 = T_BinaryOperator<
        T_TypeA1,
        T_TypeB1
    >;
    using type = CT::Vector<
        OpResult0,
        OpResult1
    >;
};

template<
    typename T_TypeA0,
    typename T_TypeA1,
    typename T_TypeA2,
    typename T_TypeB0,
    typename T_TypeB1,
    typename T_TypeB2,
    typename T_BinaryOperator
>
struct applyOperator<
    CT::Vector<
        T_TypeA0,
        T_TypeA1,
        T_TypeA2
    >,
    CT::Vector<
        T_TypeB0,
        T_TypeB1,
        T_TypeB2
    >,
    T_BinaryOperator
>
{
    using OpResult0 = T_BinaryOperator<
        T_TypeA0,
        T_TypeB0
    >;
    using OpResult1 = T_BinaryOperator<
        T_TypeA1,
        T_TypeB1
    >;
    using OpResult2 = T_BinaryOperator<
        T_TypeA2,
        T_TypeB2
    >;
    using type = CT::Vector<
        OpResult0,
        OpResult1,
        OpResult2
    >;
};

//________________________A D D____________________________

template<
    typename Lhs,
    typename Rhs
>
struct add
{
    using type = applyOperator<
        typename Lhs::vector_type,
        typename Rhs::vector_type,
        mp11::mp_plus
    >::type;
};

//________________________M U L____________________________

template<
    typename Lhs,
    typename Rhs
>
struct mul
{
    //typedef typename applyOperator<
    //typename Lhs::vector_type,
    //typename Rhs::vector_type,
    //mpl::times<mpl::_1, mpl::_2> >::type type;
    ///! Need to replace plus with times
    using type = applyOperator<
        typename Lhs::vector_type,
        typename Rhs::vector_type,
        mp11::mp_plus
    >::type;
};

//________________________M A X____________________________

template<
    typename Lhs,
    typename Rhs
>
struct max
{
    using type = applyOperator<
        typename Lhs::vector_type,
        typename Rhs::vector_type,
        mp11::mp_max
    >::type;
};

//________________________M I N____________________________

template<
    typename Lhs,
    typename Rhs = void
>
struct min
{
    using type = applyOperator<
        typename Lhs::vector_type,
        typename Rhs::vector_type,
        mp11::mp_min
    >::type;
};

template<typename T_Vec>
struct min<
    T_Vec,
    void
>
{
    using type = applyOperator<
        typename T_Vec::mplVector,
        typename T_Vec::x,
        mp11::mp_min
    >::type;
};

//________________________D O T____________________________

template<typename Lhs, typename Rhs>
struct dot
{
    /*typedef typename mul<Lhs, Rhs>::type MulResult;
    typedef typename mpl::accumulate<
        typename MulResult::mplVector,
        mp11::mp_int<0>,
        mp11::mp_plus<
            mpl::_1,
            mpl::_2
        >
    >::type type;*/
    ///! Need to replace
    using type = applyOperator<
        typename Lhs::vector_type,
        typename Rhs::vector_type,
        mp11::mp_plus
    >::type;
};

//________________________V O L U M E____________________________

//template<typename T_Vec>
//struct volume
//{
//    typedef typename mpl::accumulate<
//        typename T_Vec::mplVector,
//        mpl::int_<1>,
//        mpl::times<mpl::_1, mpl::_2>
//    >::type type;
//};

template< typename T_Vec >
struct volume
{
private:
    template< typename T_Arg >
    struct ToInt {
        static constexpr int value = T_Arg::value;
    };
    template< >
    struct ToInt< VectorParameterPlaceholder > {
        static constexpr int value = 1;
    };
public:
    using type = /*bmp11::mp_int<
        ToInt< typename At< T_Vec::mplVector, bmp11::mp_int< 0 > >::type >::value *
        ToInt< typename At< T_Vec::mplVector, bmp11::mp_int< 1 > >::type >::value *
        ToInt< typename At< T_Vec::mplVector, bmp11::mp_int< 2 > >::type >::value
    >;*/ bmp11::mp_int<10>;
};


//________________________S H R I N K T O________________________

/** shrink CT vector to given component count (dimension)
 *
 * This operation is designed to handle vectors with up to 3 components
 *
 * @tparam T_Vec vector to shrink
 * @tparam T_dim target component count
 * @treturn ::type new shrinked vector
 */
template<typename T_Vec, uint32_t T_dim>
struct shrinkTo;

template<typename T_Vec>
struct shrinkTo<T_Vec, DIM3>
{
    typedef T_Vec Vec;
    using type = CT::Vector<
        typename Vec::x,
        typename Vec::y,
        typename Vec::z
    >;
};

template<typename T_Vec>
struct shrinkTo<T_Vec, DIM2>
{
    typedef T_Vec Vec;
    using type = CT::Vector<
        typename Vec::x,
        typename Vec::y,
        VectorParameterPlaceholder
    >;
};

template<typename T_Vec>
struct shrinkTo<T_Vec, DIM1>
{
    typedef T_Vec Vec;
    using type = CT::Vector<
        typename Vec::x,
        VectorParameterPlaceholder,
        VectorParameterPlaceholder
    >;
};

//________________________A S S I G N________________________

/** Assign a type to a given component in the CT::Vector
 *
 * defines a public type as result
 *
 * @tparam T_Vec math::CT::Vector which should be changed
 * @tparam T_ComponentPos number of component to changed (type must be std::integral_constant<anyType,X>)
 * @tparam T_Value new value
 */
template<typename T_Vec, typename T_ComponentPos, typename T_Value>
struct Assign;

template<typename T_Value, typename T_0, typename T_1, typename T_2, typename T_IntegralType>
struct Assign<pmacc::math::CT::Vector<T_0, T_1, T_2>, std::integral_constant<T_IntegralType,0> , T_Value>
{
    typedef pmacc::math::CT::Vector<T_Value, T_1, T_2> type;
};

template<typename T_Value, typename T_0, typename T_1, typename T_2, typename T_IntegralType>
struct Assign<pmacc::math::CT::Vector<T_0, T_1, T_2>, std::integral_constant<T_IntegralType,1>, T_Value>
{
    typedef pmacc::math::CT::Vector<T_0, T_Value, T_2> type;
};

template<typename T_Value, typename T_0, typename T_1, typename T_2, typename T_IntegralType>
struct Assign<pmacc::math::CT::Vector<T_0, T_1, T_2>, std::integral_constant<T_IntegralType,2>, T_Value>
{
    typedef pmacc::math::CT::Vector<T_0, T_1, T_Value> type;
};

/** Assign a type to a given component in the CT::Vector if position is not out of range
 *
 * if T_ComponentPos < T_Vec::dim ? T_Value is assigned to component T_ComponentPos
 * else nothing is done.
 * defines a public type as result
 *
 * @tparam T_Vec math::CT::Vector which should be changed
 * @tparam T_ComponentPos number of component to changed (type must be std::integral_constant<anyType,X>)
 * @tparam T_Value new value
 */
template<typename T_Vec, typename T_ComponentPos, typename T_Value>
struct AssignIfInRange
{
    typedef bmp11::mp_size_t<T_Vec::dim> VectorDim;
    typedef bmp11::mp_if<
        bmp11::mp_less<T_ComponentPos, VectorDim>,
        typename pmacc::math::CT::Assign<T_Vec,T_ComponentPos,T_Value>::type,
        T_Vec
    > type;
};

//________________________At_c____________________________

/** get element from a CT::Vector
 *
 * defines a public type as result
 *
 * @tparam T_Vec input CT::Vector
 * @tparam T_idx integral index of the component
 */
template<
    typename T_Vec,
    size_t T_idx
>
struct At_c
{
    using type = mp11::mp_at_c<
        typename T_Vec::mplVector,
        T_idx
    >;
};

//________________________At____________________________

/** get element from a CT::Vector
 *
 * defines a public type as result
 *
 * @tparam T_Vec input CT::Vector
 * @tparam T_Idx integral type index of the component (e.g. boost::mp11::mp_int<2>)
 */
template<typename T_Vec, typename T_Idx>
struct At
{
    using type = mp11::mp_at<
        typename T_Vec::mplVector,
        T_Idx
    >::type;
};

//________________________make_Vector___________________

/** create CT::Vector with equal elements
 *
 * defines a public type as result
 *
 * @tparam T_dim count of components
 * @tparam T_Type type which is assigned to all components
 */
template<int T_dim, typename T_Type>
struct make_Vector;

template<typename T_Type>
struct make_Vector<1, T_Type>
{
    typedef pmacc::math::CT::Vector<T_Type> type;
};

template<typename T_Type>
struct make_Vector<2, T_Type>
{
    typedef pmacc::math::CT::Vector<T_Type, T_Type> type;
};

template<typename T_Type>
struct make_Vector<3, T_Type>
{
    typedef pmacc::math::CT::Vector<T_Type, T_Type, T_Type> type;
};

} // CT
} // math
} // pmacc
