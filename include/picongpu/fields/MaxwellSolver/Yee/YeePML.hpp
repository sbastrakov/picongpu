/* Copyright 2013-2019 Axel Huebl, Heiko Burau, Rene Widera, Benjamin Worpitz
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
#include "picongpu/fields/MaxwellSolver/Yee/Yee.hpp"
#include "picongpu/fields/MaxwellSolver/Yee/Curl.hpp"
#include "picongpu/fields/FieldManipulator.hpp"
#include "picongpu/fields/FieldE.hpp"
#include "picongpu/fields/FieldB.hpp"
#include "picongpu/fields/FieldPML.hpp"
#include "picongpu/fields/FieldManipulator.hpp"
#include "picongpu/fields/MaxwellSolver/Yee/YeePML.kernel"
#include "picongpu/fields/numericalCellTypes/NumericalCellTypes.hpp"
#include "picongpu/fields/LaserPhysics.hpp"

#include <pmacc/nvidia/functors/Assign.hpp>
#include <pmacc/mappings/threads/ThreadCollective.hpp>
#include <pmacc/memory/boxes/CachedBox.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/traits/NumberOfExchanges.hpp>

#include <cstdint>


namespace picongpu
{
namespace fields
{
namespace maxwellSolver
{

    // Yee solver with PML absorber
    // With the current design of field solvers, PML has to be specific for each solver
    // (so it is not possible to do it the same universal way the current absorber is done)
    // How to better organize it is ofc a question, maybe it is better to use trait-like design
    //
    // There probably is a way to to provide the universal PML for all solvers with a single
    // implementation, but it would require a significant redesign. Idea: elevate the current Curl
    // approach to even a higher level of abstraction, by making each solver provide the whole
    // discretization used as a type. Then a PML could be built on top of this discretization.
    // However, from the pure solver point of view this would probably look somewhat unnatural. 
    // (Probably not that relevant now.)
    template<
        typename T_CurrentInterpolation,
        class CurlE,
        class CurlB
    >
    class YeePML : public Yee< T_CurrentInterpolation, CurlE, CurlB >
    {
    private:

        using YeeSolver = Yee< T_CurrentInterpolation, CurlE, CurlB >;
        std::shared_ptr< FieldPML > fieldPML;

        // Polynomial order of the absorber strength growth towards borders
        // (often denoted 'm' or 'n' in the literature)
        uint32_t gradingOrder;

        // These parameters are specific for each side,
        // as a temporary solution sides are encoded in the absorber format (to be fixed)
        uint32_t thickness[3][2]; // unit: cells
        float3_X sigmaMax;

        void initializeParameters()
        {
            // for now reuse parameters of the old absorber
            // WARNING: this is only correct for a single process,
            // TODO: for the general case needs to take into account domain decomposition
            // so that thickness is 0 for "internal" boundaries between domains
            for (int axis = 0; axis < 3; axis++)
                for (int direction = 0; direction < 2; direction++) {
                    thickness[axis][direction] = ABSORBER_CELLS[axis][direction];
                    ///std::cout << "PML thinkess[" << axis << "][" << direction << "] = "
                    ///    << thickness[axis][direction] << std::endl;
                }

            gradingOrder = 4; // for now hardcoded with a good default value
            // Sigma optimal is based on (7.66) in Taflove 3rd ed., but uses PIC units and is divided by eps0
            float_64 sigmaOptimal[3];
            for( uint32_t dim = 0; dim < 3; dim++ )
                sigmaOptimal[ dim ] = 0.8_X * ( gradingOrder + 1 ) * SPEED_OF_LIGHT / cellSize[ dim ];
            /// todo: this should be a parameter with some default between 0.7 and 1.0
            constexpr float_X sigmaOptimalMultiplier = 1.0_X;
            for( uint32_t dim = 0; dim < 3; dim++ )
                sigmaMax[ dim ] = sigmaOptimalMultiplier * sigmaOptimal[ dim ];
        }

        bool isActiveDirection( uint32_t const exchangeIdx, uint32_t const currentStep ) const
        {
            auto const communicationMask = Environment< simDim >::get().GridController().getCommunicationMask();

            // PML can only be used on the outer borders of the simulation domain
            bool const hasNeighbor = communicationMask.isSet( exchangeIdx );
            if( hasNeighbor )
                return false;

            bool isActive = false;
            auto const relativeMask = communicationMask.getRelativeDirections< simDim >( exchangeIdx );
            for( uint32_t dim = 0; dim < simDim; dim++ )
            {
                // in relativeMask[ dim ] 1 is we are at "right" border, -1 is "left"
                int const idx = abs( relativeMask[ dim ] );
                if( thickness[ dim ][ idx ] != 0 )
                    isActive = true;
            }
            if( !isActive )
                return false;

            // this logic is copied from the old absorber,
            // after checking it seems PML should use the same logic
            /* allow to enable the absorber on the top side if the laser
            * initialization plane in y direction is *not* in cell zero
            */
            const uint32_t numSlides = MovingWindow::getInstance().getSlideCounter( currentStep );
            if( fields::laserProfiles::Selected::initPlaneY == 0 )
            {
                /* disable the absorber on top side if
                *      no slide was performed and
                *      laser init time is not over
                */
                if (numSlides == 0 && ((currentStep * DELTA_T) <= fields::laserProfiles::Selected::INIT_TIME))
                {
                    /* disable absorber on top side */
                    if( exchangeIdx == TOP )
                        return false;
                }
            }

            /* if sliding window is active we disable absorber on bottom side*/
            if (MovingWindow::getInstance().isSlidingWindowActive(currentStep) && exchangeIdx == BOTTOM)
                return false;

            return true;
        }

        void updateBHalfPML( uint32_t const currentStep )
        {
            auto const numExchanges = NumberOfExchanges< simDim >::value;
            for( uint32_t exchangeIdx = 1; exchangeIdx < numExchanges; ++exchangeIdx )
            {
                bool enabled = isActiveDirection( exchangeIdx, currentStep );
                ///std::cout << "exchangeIdx = " << exchangeIdx << ", enabled = " << enabled << "\n\n\n";
                if( !enabled )
                    continue;

                ExchangeMapping< GUARD, MappingDesc > mapper( m_cellDescription, exchangeIdx );
                constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                    pmacc::math::CT::volume< SuperCellSize >::type::value
                >::value;
                typedef SuperCellDescription<
                    SuperCellSize,
                    typename CurlE::LowerMargin,
                    typename CurlE::UpperMargin
                > BlockArea;

                yeePML::detail::Parameters parameters;
                parameters.normalizedSigmaMax = this->sigmaMax;
                parameters.gradingOrder = this->gradingOrder;
                PMACC_KERNEL( yeePML::KernelUpdateBHalf< numWorkers, BlockArea >{ } )
                    ( mapper.getGridDim(), numWorkers )(
                        CurlE(),
                        this->fieldPML->getDeviceDataBox(),
                        this->fieldB->getDeviceDataBox(),
                        this->fieldE->getDeviceDataBox(),
                        mapper,
                        parameters
                    );
            }
        }

        void updateEPML( uint32_t currentStep )
        {
            auto const numExchanges = NumberOfExchanges< simDim >::value;
            for( uint32_t exchangeIdx = 1; exchangeIdx < numExchanges; ++exchangeIdx )
            {
                bool enabled = isActiveDirection( exchangeIdx, currentStep );
                ///std::cout << "exchangeIdx = " << exchangeIdx << ", enabled = " << enabled << "\n\n\n";
                if( !enabled )
                    continue;

                ExchangeMapping< GUARD, MappingDesc > mapper( m_cellDescription, exchangeIdx );
                constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                    pmacc::math::CT::volume< SuperCellSize >::type::value
                >::value;
                typedef SuperCellDescription<
                    SuperCellSize,
                    typename CurlB::LowerMargin,
                    typename CurlB::UpperMargin
                > BlockArea;

                yeePML::detail::Parameters parameters;
                parameters.normalizedSigmaMax = this->sigmaMax;
                parameters.gradingOrder = this->gradingOrder;
                PMACC_KERNEL( yeePML::KernelUpdateE< numWorkers, BlockArea >{ } )
                    ( mapper.getGridDim(), numWorkers )(
                        CurlE(),
                        this->fieldPML->getDeviceDataBox(),
                        this->fieldE->getDeviceDataBox(),
                        this->fieldB->getDeviceDataBox(),
                        mapper,
                        parameters
                    );
            }
        }

    public:

        YeePML(MappingDesc cellDescription) : Yee(cellDescription)
        {
            DataConnector &dc = Environment<>::get().DataConnector();
            this->fieldPML = dc.get< FieldPML >( FieldPML::getName(), true );
            initializeParameters();
        }

        // TODO: the following update_beforeCurrent() and update_afterCurrent()
        // are exactly the same as in the Yee solver (but internally called
        // updateBHalf, updateE are not). It would be ofc much better to not
        // duplicate this logic, and this should be possible rather easily
        void update_beforeCurrent( uint32_t const currentStep )
        {
            // here we need up-to-date E everywhere, including PML for near-PML border
            YeeSolver::updateBHalf< CORE + BORDER >();
            updateBHalfPML( currentStep );
            EventTask eRfieldB = fieldB->asyncCommunication( __getTransactionEvent() );

            YeeSolver::updateE< CORE >();
            __setTransactionEvent( eRfieldB );
            YeeSolver::updateE< BORDER >();
        }

        void update_afterCurrent( uint32_t const currentStep )
        {
            updateEPML( currentStep );
            if( laserProfiles::Selected::INIT_TIME > 0.0_X )
                LaserPhysics{}( currentStep );

            EventTask eRfieldE = fieldE->asyncCommunication( __getTransactionEvent() );

            YeeSolver::updateBHalf< CORE >();
            __setTransactionEvent( eRfieldE );
            YeeSolver::updateBHalf< BORDER >();
            updateBHalfPML( currentStep );

            EventTask eRfieldB = fieldB->asyncCommunication( __getTransactionEvent() );
            __setTransactionEvent( eRfieldB );
        }

        static pmacc::traits::StringProperty getStringProperties()
        {
            pmacc::traits::StringProperty propList( "name", "YeePML" );
            return propList;
        }
    };

} // namespace maxwellSolver
} // namespace fields
} // picongpu
