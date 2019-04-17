/* Copyright 2013-2019 Axel Huebl, Heiko Burau, Rene Widera, Felix Schmitt,
 *                     Richard Pausch, Benjamin Worpitz
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

#include "picongpu/fields/MaxwellSolver/YeePML/SplitFields.hpp"

#include <pmacc/eventSystem/EventSystem.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>
#include <pmacc/mappings/kernel/ExchangeMapping.hpp>
#include <pmacc/memory/buffers/GridBuffer.hpp>

#include "picongpu/fields/FieldManipulator.hpp"

#include <pmacc/dimensions/SuperCellDescription.hpp>

#include "picongpu/fields/MaxwellSolver/Solvers.hpp"
#include "picongpu/fields/numericalCellTypes/NumericalCellTypes.hpp"

#include <pmacc/math/Vector.hpp>

#include "picongpu/particles/traits/GetInterpolation.hpp"
#include <pmacc/particles/traits/FilterByFlag.hpp>

#include "picongpu/traits/GetMargin.hpp"
#include "picongpu/traits/SIBaseUnits.hpp"
#include "picongpu/particles/traits/GetMarginPusher.hpp"

#include <boost/mpl/accumulate.hpp>

#include <list>
#include <iostream>
#include <memory>


namespace picongpu
{
namespace fields
{
namespace maxwellSolver
{
namespace yeePML
{

SplitFields::SplitFields( MappingDesc cellDescription ) :
SimulationFieldHelper<MappingDesc>( cellDescription )
{
    data.reset(
        new GridBuffer<ValueType, simDim > ( cellDescription.getGridLayout( ) )
    );
}

SimulationDataId SplitFields::getUniqueId()
{
    return getName();
}

void SplitFields::synchronize( )
{
    data->deviceToHost( );
}

void SplitFields::syncToDevice( )
{
    data->hostToDevice( );
}

EventTask SplitFields::asyncCommunication( EventTask serialEvent )
{
    EventTask eB = data->asyncCommunication( serialEvent );
    return eB;
}

GridLayout<simDim> SplitFields::getGridLayout( )
{

    return cellDescription.getGridLayout( );
}

SplitFields::DataBoxType SplitFields::getHostDataBox( )
{

    return data->getHostBuffer( ).getDataBox( );
}

SplitFields::DataBoxType SplitFields::getDeviceDataBox( )
{

    return data->getDeviceBuffer( ).getDataBox( );
}

GridBuffer<SplitFields::ValueType, simDim> &SplitFields::getGridBuffer( )
{

    return *data;
}

void SplitFields::reset( uint32_t )
{
    data->getHostBuffer( ).reset( true );
    data->getDeviceBuffer( ).reset( false );
}

HDINLINE
SplitFields::UnitValueType
SplitFields::getUnit( )
{
    return UnitValueType( 1.0_X );
}

HINLINE
std::vector<float_64>
SplitFields::getUnitDimension( )
{
    std::vector<float_64> unitDimension( 7, 0.0 );
    return unitDimension;
}

std::string
SplitFields::getName( )
{
    return "PML split fields";
}

uint32_t
SplitFields::getCommTag( )
{
    // These fields do not need to be communicated
    return 0;
}

} // namespace yeePML
} // namespace maxwellSolver
} // namespace fields
} // namespace picongpu
