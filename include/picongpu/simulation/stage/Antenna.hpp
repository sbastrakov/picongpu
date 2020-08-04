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

#include "picongpu/fields/antenna/ApplyAntenna.hpp"
#include "picongpu/fields/FieldJ.hpp"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/Environment.hpp>

#include <cstdint>


namespace picongpu
{
namespace simulation
{
namespace stage
{

    //! Functor for the stage of the PIC loop applying antenna to generate current
    class Antenna
    {
    public:

        /** Create an antenna functor
         *
         * Having this in constructor is a temporary solution.
         *
         * @param cellDescription mapping for kernels
         */
        Antenna( MappingDesc const cellDescription ):
            cellDescription( cellDescription )
        {
        }

        /** Add the antenna current to the current density
         *
         * @param step index of time iteration
         */
        void operator( )( uint32_t const step ) const
        {
            using namespace pmacc;
            DataConnector & dc = Environment< >::get( ).DataConnector( );
            auto & fieldJ = *dc.get< FieldJ >( FieldJ::getName( ), true );
            fields::antenna::ApplyAntenna< fields::Antenna > applyAntenna;
            applyAntenna(
                fieldJ,
                step,
                cellDescription
            );
            dc.releaseData( FieldJ::getName( ) );
        }

    private:

        //! Mapping for kernels
        MappingDesc cellDescription;

    };

} // namespace stage
} // namespace simulation
} // namespace picongpu
