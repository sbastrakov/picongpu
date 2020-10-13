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

#include "picongpu/fields/incidentField/Profiles.def"
#include "picongpu/simulation_defines.hpp"


namespace picongpu
{
namespace fields
{
namespace incidentField
{

    //! Concept defining interface of incident field functors for E and B
    struct FunctorIncidentFieldConcept
    {
        /** Create a functor
         *
         * @param unitField conversion factor from SI to internal units,
         *                  field_internal = field_SI / unitField
         */
        HDINLINE FunctorIncidentFieldConcept( float3_64 unitField );

        /** Return incident field for the given position and time.
         *
         * Note that component along the source boundary (e.g. y component for
         * YMin and YMax boundaries) will not be used, and so should be zero.
         *
         * @param cellIdx cell index in global domain, note that it is fractional
         * @param currentStep current time step index, note that it is fractional
         * @return incident field value in internal units
         */
        HDINLINE float3_X operator()(
            floatD_X const & cellIdx,
            float_X currentStep
        ) const;
    };

    template<
        typename T_FunctorIncidentE,
        typename T_FunctorIncidentB
    >
    struct Source
    {
        //! Incident E functor
        using FunctorIncidentE = T_FunctorIncidentE;

        //! Incident B functor
        using FunctorIncidentB = T_FunctorIncidentB;
    };

} // namespace incidentField
} // namespace fields
} // namespace picongpu
