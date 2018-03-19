/* Copyright 2013-2018 Axel Huebl, Heiko Burau, Rene Widera, Benjamin Worpitz
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
#include "pmacc/forward.hpp"

#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/utility.hpp>


namespace pmacc
{
namespace algorithms
{
namespace forEach
{
namespace detail
{
    /** call the functor were itBegin points to
     *
     *  \tparam itBegin iterator to an element in a mpl sequence
     *  \tparam itEnd iterator to the end of a mpl sequence
     *  \tparam isEnd true if itBegin == itEnd, else false
     */
    template<
        typename List,
        bool isEnd = bmp11::is_empty< List >::value
    >
    struct CallFunctorOfIterator
    {
        PMACC_NO_NVCC_HDWARNING
        template< typename ... T_Types >
        HDINLINE void
        operator( )( T_Types const & ... ts ) const
        {
            bmp11::mp_front< List >( )( getForwardedValue( ts ) ... );
            bmp11::mp_pop_front< List >( )( ts ... );
        }

    };

    /** Recursion end of ForEach */
    template< typename List >
    struct CallFunctorOfIterator<
        List,
        true
    >
    {
        PMACC_NO_NVCC_HDWARNING
        template< typename ... T_Types >
        HDINLINE void
        operator()( T_Types const & ... ) const
        {
        }
    };

} // namespace detail

    /** Compile-Time for each for Boost::MPL Type Lists
     *
     *  \tparam T_MPLSeq A mpl sequence that can be accessed by mpl::begin, mpl::end, mpl::next
     *  \tparam T_Functor An unary lambda functor with a HDINLINE void operator()(...) method
     *          _1 is substituted by Accessor's result.
     *          The maximum number of parameters for the operator() is limited by
     *          PMACC_MAX_FUNCTOR_OPERATOR_PARAMS
     *  \tparam T_Accessor An unary lambda operation
     *
     * Example:
     *      MPLSeq = boost::mp11::mp_list<int,float>
     *      Functor = any unary lambda functor
     *      Accessor = lambda operation identity
     *
     *      definition: F(X) means F<X>
     *
     *      call:   ForEach<MPLSeq,Functor,Accessor>()(42);
     *      unrolled code: Functor(Accessor(int))(42);
     *                     Functor(Accessor(float))(42);
     */
    template<
        typename T_MPLSeq,
        typename T_Functor,
        typename T_Accessor = compileTime::accessors::Identity< >
    >
    struct ForEach
    {

        template< typename X >
        struct ReplacePlaceholder
        {
            using type = T_Functor< typename T_Accessor< X >::type >::type;
        };

        using SolvedFunctors = bmp11::mp_transform<
            bmp11::mp_identity_t< ReplacePlaceholder >,
            T_MPLSeq
        >;

        PMACC_NO_NVCC_HDWARNING
        template< typename ... T_Types >
        HDINLINE void
        operator( )( const T_Types& ... ts ) const
        {
            detail::CallFunctorOfIterator< SolvedFunctors >::NextCall()( ts ... );
        }

    };

} // namespace forEach
} // namespace algorithms
} // namespace pmacc
