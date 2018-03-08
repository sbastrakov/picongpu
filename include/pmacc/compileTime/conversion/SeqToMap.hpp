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

#include "pmacc/types.hpp"
#include "pmacc/compileTime/accessors/Identity.hpp"

#include <boost/mpl/map.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/bind.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/type_traits.hpp>


namespace pmacc
{

/** convert boost mpl sequence to a mpl map
 *
 * @tparam T_MPLSeq any boost mpl sequence
 * @tparam T_UnaryOperator unary operator to translate type from the sequence
 * to a mpl pair
 * @tparam T_Accessor An unary lambda operator which is used before the type
 * from the sequence is passed to T_UnaryOperator
 * @return ::type mpl map
 */
template<typename T_MPLSeq,
typename T_UnaryOperator,
typename T_Accessor = compileTime::accessors::Identity<>
>
struct SeqToMap
{

    template<typename X>
    struct Op : typename T_UnaryOperator< typename T_Accessor< X >::type >::type
    {
    };

    typedef T_MPLSeq MPLSeq;
    typedef bmpl::inserter< bmpl::map<>, bmpl::insert<bmp11::_1, bmp11::_2> > Map_inserter;
    using type = bmp11::mp_transform<
        MPLSeq,
        Op<bmp11::_1> ,
        Map_inserter
    >;
};

}//namespace pmacc
