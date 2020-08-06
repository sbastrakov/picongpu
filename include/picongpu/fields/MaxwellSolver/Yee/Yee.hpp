/* Copyright 2013-2020 Axel Huebl, Heiko Burau, Rene Widera, Benjamin Worpitz
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
#include "picongpu/fields/MaxwellSolver/Yee/Yee.def"
#include "picongpu/fields/absorber/ExponentialDamping.hpp"
#include "picongpu/fields/antenna/ApplyAntenna.hpp"
#include "picongpu/fields/FieldE.hpp"
#include "picongpu/fields/FieldB.hpp"
#include "picongpu/fields/MaxwellSolver/Yee/Yee.kernel"
#include "picongpu/fields/cellType/Yee.hpp"
#include "picongpu/fields/LaserPhysics.hpp"
#include "picongpu/traits/GetMargin.hpp"

#include <pmacc/nvidia/functors/Assign.hpp>
#include <pmacc/mappings/threads/ThreadCollective.hpp>
#include <pmacc/memory/boxes/CachedBox.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>


namespace picongpu
{
namespace fields
{
namespace maxwellSolver
{

    template<
        typename T_CurrentInterpolation,
        class CurlE,
        class CurlB
    >
    class Yee
    {
    private:
        typedef MappingDesc::SuperCellSize SuperCellSize;


        std::shared_ptr< FieldE > fieldE;
        std::shared_ptr< FieldB > fieldB;
        MappingDesc m_cellDescription;

        template<uint32_t AREA>
        void updateE()
        {
            /* Courant-Friedrichs-Levy-Condition for Yee Field Solver:
             *
             * A workaround is to add a template dependency to the expression.
             * `sizeof(ANY_TYPE*) != 0` is always true and defers the evaluation.
             */
            PMACC_CASSERT_MSG(Courant_Friedrichs_Levy_condition_failure____check_your_grid_param_file,
                (SPEED_OF_LIGHT*SPEED_OF_LIGHT*DELTA_T*DELTA_T*INV_CELL2_SUM)<=1.0 && sizeof(T_CurrentInterpolation*) != 0);

            typedef SuperCellDescription<
                    SuperCellSize,
                    typename traits::GetLowerMargin<CurlB>::type,
                    typename traits::GetUpperMargin<CurlB>::type
                    > BlockArea;

            AreaMapping<AREA, MappingDesc> mapper(m_cellDescription);

            constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                pmacc::math::CT::volume< SuperCellSize >::type::value
            >::value;

            PMACC_KERNEL(yee::KernelUpdateE< numWorkers, BlockArea >{ })
                ( mapper.getGridDim(), numWorkers )(
                    CurlB( ),
                    this->fieldE->getDeviceDataBox(),
                    this->fieldB->getDeviceDataBox(),
                    mapper
                );
        }

        template<uint32_t AREA>
        void updateBHalf()
        {
            typedef SuperCellDescription<
                    SuperCellSize,
                    typename CurlE::LowerMargin,
                    typename CurlE::UpperMargin
                    > BlockArea;

            AreaMapping<AREA, MappingDesc> mapper(m_cellDescription);

            constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                pmacc::math::CT::volume< SuperCellSize >::type::value
            >::value;

            PMACC_KERNEL(yee::KernelUpdateBHalf< numWorkers, BlockArea >{ })
                ( mapper.getGridDim(), numWorkers )(
                    CurlE( ),
                    this->fieldB->getDeviceDataBox(),
                    this->fieldE->getDeviceDataBox(),
                    mapper
                );
        }

    public:

        using CellType = cellType::Yee;
        using CurrentInterpolation = T_CurrentInterpolation;

        Yee(MappingDesc cellDescription) : m_cellDescription(cellDescription)
        {
            DataConnector &dc = Environment<>::get().DataConnector();

            this->fieldE = dc.get< FieldE >( FieldE::getName(), true );
            this->fieldB = dc.get< FieldB >( FieldB::getName(), true );
        }

        void update_beforeCurrent(uint32_t currentStep)
        {
            updateBHalf < CORE+BORDER >();
            fields::antenna::ApplyAntenna applyAntenna;
            // update B by half step, to step = currentStep + 0.5, so step for E_inc = currentStep
            applyAntenna.updateB(
                static_cast< float_X >( currentStep ),
                0.5_X,
                m_cellDescription
            );
            EventTask eRfieldB = fieldB->asyncCommunication(__getTransactionEvent());

            updateE<CORE>();
            // Antenna update does not use exchanged B, so does not have to wait for it
            // update E by full step, to step = currentStep + 1, so step for B_inc = currentStep + 0.5
            applyAntenna.updateE(
                static_cast< float_X >( currentStep ) + 0.5_X,
                1.0_X,
                m_cellDescription
            );
            __setTransactionEvent(eRfieldB);
            updateE<BORDER>();
        }

        void update_afterCurrent(uint32_t currentStep)
        {
            using Absorber = absorber::ExponentialDamping;
            Absorber::run(
                currentStep,
                this->m_cellDescription,
                this->fieldE->getDeviceDataBox()
            );
            if (laserProfiles::Selected::INIT_TIME > float_X(0.0))
                LaserPhysics{}(currentStep);

            // Antenna update does not use exchanged E, so does not have to wait for it
            fields::antenna::ApplyAntenna applyAntenna;
            // update B by half step, to step currentStep + 1.5, so step for E_inc = currentStep + 1
            applyAntenna.updateB(
                static_cast< float_X >( currentStep ) + 1.0_X,
                0.5_X,
                m_cellDescription
            );

            EventTask eRfieldE = fieldE->asyncCommunication(__getTransactionEvent());

            updateBHalf < CORE> ();
            __setTransactionEvent(eRfieldE);
            updateBHalf < BORDER > ();

            Absorber::run(
                currentStep,
                this->m_cellDescription,
                fieldB->getDeviceDataBox()
            );

            EventTask eRfieldB = fieldB->asyncCommunication(__getTransactionEvent());
            __setTransactionEvent(eRfieldB);
        }

        static pmacc::traits::StringProperty getStringProperties()
        {
            pmacc::traits::StringProperty propList( "name", "Yee" );
            return propList;
        }
    };

} // namespace maxwellSolver
} // namespace fields
} // picongpu
