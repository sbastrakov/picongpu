/* Copyright 2014-2018 Rene Widera
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

#include "pmacc/HandleGuardRegion.hpp"
#include "pmacc/particles/policies/ExchangeParticles.hpp"
#include "pmacc/particles/policies/DeleteParticles.hpp"
#include "pmacc/compileTime/conversion/ToSeq.hpp"

#include <boost/mp11/list.hpp>


namespace pmacc
{

/** ParticleDescription defines attributes, methods and flags of a particle
 *
 * This class holds no runtime data.
 * The class holds information about the name, attributes, flags and methods of a
 * particle.
 *
 * @tparam T_Name name of described particle (e.g. electron, ion)
 *                type must be a boost::mpl::string
 * @tparam T_SuperCellSize compile time size of a super cell
 * @tparam T_ValueTypeSeq sequence or single type with value_identifier
 * @tparam T_Flags sequence or single type with identifier to add flags on a frame
 * @tparam T_MethodsList sequence or single class with particle methods
 *                       (e.g. calculate mass, gamma, ...)
 *                       (e.g. useSolverXY, calcRadiation, ...)
 * @tparam T_FrameExtensionList sequence or single class with frame extensions
 *                    - extension must be an unary template class that defines type
 *                    - type of the final frame is applied to each extension class
 *                      (this allows pointers and references to a frame itself)
 *                    - the final frame that uses ParticleDescription inherits from all
 *                      extension classes
 */
template<
typename T_Name,
typename T_SuperCellSize,
typename T_ValueTypeSeq,
typename T_Flags = bmp11::mp_list<>,
typename T_HandleGuardRegion = HandleGuardRegion<particles::policies::ExchangeParticles, particles::policies::DeleteParticles>,
typename T_MethodsList = bmp11::mp_list<>,
typename T_FrameExtensionList = bmp11::mp_list<>
>
struct ParticleDescription
{
    typedef T_Name Name;
    typedef T_SuperCellSize SuperCellSize;
    using ValueTypeSeq = ToSeq<T_ValueTypeSeq>;
    using FlagsList = ToSeq<T_Flags>;
    typedef T_HandleGuardRegion HandleGuardRegion;
    using MethodsList = ToSeq<T_MethodsList>;
    using FrameExtensionList = ToSeq<T_FrameExtensionList>;
    typedef ParticleDescription<
        Name,
        SuperCellSize,
        ValueTypeSeq,
        FlagsList,
        HandleGuardRegion,
        MethodsList,
        FrameExtensionList
    > ThisType;

};


/** Get ParticleDescription with a new ValueTypeSeq
 *
 * @tparam T_OldParticleDescription base description
 * @tparam T_NewValueTypeSeq new boost mpl sequence with value types
 * @treturn ::type new ParticleDescription
 */
template<typename T_OldParticleDescription, typename T_NewValueTypeSeq>
struct ReplaceValueTypeSeq
{
    typedef T_OldParticleDescription OldParticleDescription;
    typedef ParticleDescription<
    typename OldParticleDescription::Name,
    typename OldParticleDescription::SuperCellSize,
    ToSeq<T_NewValueTypeSeq>,
    typename OldParticleDescription::FlagsList,
    typename OldParticleDescription::HandleGuardRegion,
    typename OldParticleDescription::MethodsList,
    typename OldParticleDescription::FrameExtensionList
    > type;
};

/** Get ParticleDescription with a new FrameExtensionSeq
 *
 * @tparam T_OldParticleDescription base description
 * @tparam T_FrameExtensionSeq new boost mpl sequence with value types
 * @treturn ::type new ParticleDescription
 */
template<typename T_OldParticleDescription, typename T_FrameExtensionSeq>
struct ReplaceFrameExtensionSeq
{
    typedef T_OldParticleDescription OldParticleDescription;
    typedef ParticleDescription<
    typename OldParticleDescription::Name,
    typename OldParticleDescription::SuperCellSize,
    typename OldParticleDescription::ValueTypeSeq,
    typename OldParticleDescription::FlagsList,
    typename OldParticleDescription::HandleGuardRegion,
    typename OldParticleDescription::MethodsList,
    ToSeq<T_FrameExtensionSeq>
    > type;
};

} //namespace pmacc
