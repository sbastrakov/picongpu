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


#include <boost/mp11/algorithm.hpp>


namespace pmacc
{

/* remove types from a list
 *
 * @tparam T_SeqSrc source mp11 list from were we delete types
 * @tparam T_SeqObjectsToRemove mp11 list with types which should be deleted
 */
template<
    typename T_SeqSrc,
    typename T_SeqObjectsToRemove
>
struct RemoveFromSeq
{
private:
    template< typename T_Value >
    using HasId = bmp11::mp_contains<
        T_SeqObjectsToRemove,
        T_Value
    >;

public:
    using type = bmp11::mp_remove_if<
        T_SeqSrc,
        HasId
    >;
};

}//namespace pmacc
