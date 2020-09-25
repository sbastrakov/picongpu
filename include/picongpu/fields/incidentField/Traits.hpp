/* Copyright 2020 Sergei Bastrakov
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "picongpu/simulation_defines.hpp"

#include <picongpu/fields/incidentField/Profiles.hpp>

#include <cstdint>


namespace picongpu
{
namespace fields
{
namespace incidentField
{
namespace detail
{

    /** Get ZMin incident field source type depending on the 2d/3d case
     *
     * The resulting field source is set as ::type.
     * Implementation for 3d returns ZMin from incidentField.param.
     *
     * @tparam T_is3d is the simulation 3d
     */
    template< bool T_is3d = true >
    struct GetZMin
    {
        using type = ZMin;
    };

    //! Get None incident field source type as ZMin in the 2d case
    template< >
    struct GetZMin< false >
    {
        using type = None;
    };

    /** Get ZMax incident field source type depending on the 2d/3d case
     *
     * The resulting field source is set as ::type.
     * Implementation for 3d returns ZMin from incidentField.param.
     *
     * @tparam T_is3d is the simulation 3d
     */
    template< bool T_is3d = true >
    struct GetZMax
    {
        using type = ZMax;
    };

    //! Get None incident field source type as ZMax in the 2d case
    template< >
    struct GetZMax< false >
    {
        using type = None;
    };

} // namespace detail

    /** Get ZMin incident field source type
     *
     * The resulting field source is set as ::type.
     * In 3d returns ZMin from incidentField.param, in 2d returns None.
     */
    struct GetZMin
    {
        using type = detail::GetZMin< simDim == 3 >::type;
    };

    /** Get ZMax incident field source type
     *
     * The resulting field source is set as ::type.
     * In 3d returns ZMax from incidentField.param, in 2d returns None.
     */
    struct GetZMax
    {
        using type = detail::GetZMax< simDim == 3 >::type;
    };

} // namespace incidentField
} // namespace fields
} // namespace picongpu
