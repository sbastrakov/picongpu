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

#include "picongpu/fields/FieldPML.hpp"

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

using namespace pmacc;

FieldPML::FieldPML( MappingDesc cellDescription ) :
SimulationFieldHelper<MappingDesc>( cellDescription )
{
    /*#####create FieldB###############*/
    fieldPML = new GridBuffer<ValueType, simDim > ( cellDescription.getGridLayout( ) );
}

FieldPML::~FieldPML( )
{
    __delete(fieldB);
}

SimulationDataId FieldPML::getUniqueId()
{
    return getName();
}

void FieldPML::synchronize( )
{
    fieldPML->deviceToHost( );
}

void FieldPML::syncToDevice( )
{

    fieldPML->hostToDevice( );
}

EventTask FieldPML::asyncCommunication( EventTask serialEvent )
{

    EventTask eB = fieldPML->asyncCommunication( serialEvent );
    return eB;
}

GridLayout<simDim> FieldPML::getGridLayout( )
{

    return cellDescription.getGridLayout( );
}

FieldPML::DataBoxType FieldPML::getHostDataBox( )
{

    return fieldPML->getHostBuffer( ).getDataBox( );
}

FieldPML::DataBoxType FieldPML::getDeviceDataBox( )
{

    return fieldPML->getDeviceBuffer( ).getDataBox( );
}

GridBuffer<FieldPML::ValueType, simDim> &FieldPML::getGridBuffer( )
{

    return *fieldPML;
}

void FieldPML::reset( uint32_t )
{
    fieldPML->getHostBuffer( ).reset( true );
    fieldPML->getDeviceBuffer( ).reset( false );
}

HDINLINE
FieldPML::UnitValueType
FieldPML::getUnit( )
{
    return UnitValueType( 1.0_X );
}

HINLINE
std::vector<float_64>
FieldPML::getUnitDimension( )
{
    std::vector<float_64> unitDimension( 7, 0.0 );
    return unitDimension;
}

std::string
FieldPML::getName( )
{
    return "PML";
}

uint32_t
FieldPML::getCommTag( )
{
    return 0; /// todo
}

} //namespace picongpu
