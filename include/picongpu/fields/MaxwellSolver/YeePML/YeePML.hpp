/* Copyright 2019 Sergei Bastrakov
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
#include "picongpu/fields/MaxwellSolver/YeePML/Solver.hpp"
#include "picongpu/fields/numericalCellTypes/NumericalCellTypes.hpp"

#include <pmacc/traits/GetStringProperties.hpp>


namespace picongpu
{
namespace fields
{
namespace maxwellSolver
{

    /**
     * Yee solver with split-field perfectly matched layer (PML).
     * It implements the field solver interface and is supposed to be
     * used from param files and in the computatioanl loop.
     *
     * It supports use model of new instance on every time step.
     * However, PML solver has a state that has to be preserved between time
     * iterations. Therefore, this solver is a wrapper around a singleton of
     * yeePML::Solver.
     */
    template<
        typename T_CurrentInterpolation,
        typename T_CurlE,
        typename T_CurlB
    >
    class YeePML
    {
    public:

        using NummericalCellType = picongpu::numericalCellTypes::YeeCell;
        using CurrentInterpolation = T_CurrentInterpolation;

        YeePML(MappingDesc cellDescription)
        {
            // Initialize the solver at first call
            auto& solverPtr = getSolverPtr();
            if( !solverPtr )
                solverPtr.reset( new Solver( cellDescription ) );
        }

        void update_beforeCurrent( uint32_t const currentStep )
        {
            getSolverPtr()->update_beforeCurrent( currentStep );
        }

        void update_afterCurrent( uint32_t const currentStep )
        {
            getSolverPtr()->update_afterCurrent( currentStep );
        }

        static pmacc::traits::StringProperty getStringProperties()
        {
            pmacc::traits::StringProperty propList( "name", "YeePML" );
            return propList;
        }

    private:

        using Solver = yeePML::Solver<
            T_CurrentInterpolation,
            T_CurlE,
            T_CurlB
        >;

        static std::unique_ptr<Solver>& getSolverPtr()
        {
            static std::unique_ptr<Solver> instance;
            return instance;
        }

    };

} // namespace maxwellSolver
} // namespace fields
} // namespace picongpu
