/* Copyright 2019 Axel Huebl, Heiko Burau, Rene Widera, Benjamin Worpitz,
 *                Sergei Bastrakov
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
#include "picongpu/fields/MaxwellSolver/YeePML/YeePML.kernel"
#include "picongpu/fields/numericalCellTypes/NumericalCellTypes.hpp"

#include <pmacc/traits/GetStringProperties.hpp>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>


namespace picongpu
{
namespace fields
{
namespace maxwellSolver
{

    /**
      * Yee solver with perfectly matched layer (PML) absorber.
      */
    template<
        typename T_CurrentInterpolation,
        typename T_CurlE,
        typename T_CurlB
    >
    class YeePML
    {
    private:

        // Yee solver data
        std::shared_ptr< FieldE > fieldE;
        std::shared_ptr< FieldB > fieldB;
        MappingDesc cellDescription;

        // PML data
        std::shared_ptr< yeePML::SplitFields > splitFields;
        /** Thickness of the layer in number of cells
         */
        struct Thickness
        {
            DataSpace< simDim > negativeBorder;
            DataSpace< simDim > positiveBorder;

            /** Provides element access with indexing of the old absorber:
             *  axis is 0 = x, 1 = y, or 2 = z; direction is 0 = negative, 1 = positive
             *  This is only for initialization convenience and so does not have a device version
             */
            int& operator()( uint32_t const axis, uint32_t const direction )
            {
                // Since this is not performance-critical at all, do range checks
                if( axis >= simDim )
                    throw std::out_of_range("In Thickness::operator() the axis = " +
                        std::to_string( axis ) + " is invalid");
                if( direction == 0 )
                    return negativeBorder[ axis ];
                else
                    if( direction == 1 )
                        return positiveBorder[ axis ];
                    else
                        throw std::out_of_range("In Thickness::operator() the direction = " +
                            std::to_string( direction ) +  " is invalid");
            }
        };

        /* Thickness in terms of the global domain
         * We store only global thickness, as the local one can change
         * during the simulation
         */
        Thickness globalSize;

        // Strength-related parameters
        float3_X sigmaMax;
        // Polynomial order of the absorber strength growth towards borders
        // (often denoted 'm' or 'n' in the literature)
        uint32_t gradingOrder;

        void initializeParameters( )
        {
            globalSize = getGlobalThickness( );
            gradingOrder = 4; // for now hardcoded with a good default value
            // Sigma optimal is based on (7.66) in Taflove 3rd ed.,
            // but uses PIC units and is divided by eps0
            float_64 sigmaOptimal[3];
            for( auto dim = 0; dim < 3; dim++ )
                sigmaOptimal[ dim ] = 0.8_X * ( gradingOrder + 1 ) * SPEED_OF_LIGHT / cellSize[ dim ];
            // This can become a a parameter with some default between 0.7 and 1.0
            constexpr float_X sigmaOptimalMultiplier = 1.0_X;
            for( auto dim = 0; dim < 3; dim++ )
                sigmaMax[ dim ] = sigmaOptimalMultiplier * sigmaOptimal[ dim ];
        }

        Thickness getGlobalThickness( ) const
        {
            Thickness globalThickness;
            // For now thickness of the exponential damping absorber is used
            for( auto axis = 0; axis < simDim; axis++ )
                for( auto direction = 0; direction < 2; direction++ )
                    globalThickness( axis, direction ) = ABSORBER_CELLS[ axis ][ direction ];
            return globalThickness;
        }

        /** Get local PML thickness for the current time step
         *  It depends on the current step because of the moving window
         */
        Thickness getLocalThickness( uint32_t const currentStep ) const
        {
            auto & movingWindow = MovingWindow::getInstance();
            auto const numSlides = movingWindow.getSlideCounter( currentStep );
            auto const numExchanges = NumberOfExchanges< simDim >::value;
            auto const communicationMask = Environment< simDim >::get().GridController().getCommunicationMask();
            // The logic of these checks is the same as in FieldManipulator::absorbBorder(),
            // to disable the absorber at a border we set the corresponding thickness to 0
            Thickness localThickness = globalSize;
            for( auto exchange = 1; exchange < numExchanges; ++exchange )
            {
                // Ignore directions except left, right, top, bottom, back, front
                if( FRONT % exchange != 0 )
                    continue;

                // Transform exchange into a pair of axis and direction
                uint32_t axis = 0;
                if( exchange >= BOTTOM && exchange <= TOP )
                    axis = 1;
                if( exchange >= BACK )
                    axis = 2;
                uint32_t direction = exchange % 2;

                // No PML at the borders between two local domains
                bool hasNeighbour = communicationMask.isSet( exchange );
                if( hasNeighbour )
                    localThickness( axis, direction ) = 0;

                // Disable PML during laser initialization
                if( fields::laserProfiles::Selected::initPlaneY == 0 )
                {
                    bool isLaserInitializationOver =
                        (currentStep * DELTA_T) > fields::laserProfiles::Selected::INIT_TIME;
                    if( numSlides == 0 && !isLaserInitializationOver && exchange == TOP )
                        localThickness( axis, direction ) = 0;
                }

                // Disable PML at the far side of the moving window
                if( movingWindow.isSlidingWindowActive( currentStep ) && exchange == BOTTOM )
                    localThickness( axis, direction ) = 0;
            }
            return localThickness;
        }

        yeePML::detail::Parameters computeParameters( uint32_t const currentStep ) const
        {
            Thickness localThickness = getLocalThickness( currentStep );
            yeePML::detail::Parameters parameters;
            // Convert size type here to avoid doing that in kernels
            for( auto axis = 0; axis < simDim; axis++ )
            {
                parameters.negativeBorderSize[ axis ] = static_cast< float_X >(
                    localThickness.negativeBorder[ axis ]
                );
                parameters.positiveBorderSize[ axis ] = static_cast< float_X >(
                    localThickness.positiveBorder[ axis ]
                );
            }
            parameters.normalizedSigmaMax = this->sigmaMax;
            parameters.gradingOrder = this->gradingOrder;
            return parameters;
        }

        template< uint32_t AREA >
        void updateBHalf( uint32_t const currentStep )
        {
            typedef SuperCellDescription<
                SuperCellSize,
                typename CurlE::LowerMargin,
                typename CurlE::UpperMargin
            > BlockArea;
            AreaMapping< AREA, MappingDesc > mapper( cellDescription );
            constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                pmacc::math::CT::volume< SuperCellSize >::type::value
            >::value;
            // Note: here it is possible to check if PML is enabled in the local
            // domain at all, and otherwise optimize by calling the normal Yee update,
            // but we do not do that, as this is fragile wrt plans for close future changes
            PMACC_KERNEL(yeePML::KernelUpdateBHalf< numWorkers, BlockArea >{ })
                ( mapper.getGridDim(), numWorkers )(
                    CurlE( ),
                    this->splitFields->getDeviceDataBox(),
                    this->fieldB->getDeviceDataBox(),
                    this->fieldE->getDeviceDataBox(),
                    mapper,
                    computeParameters( currentStep )
                );
        }

        template< uint32_t AREA >
        void updateE( uint32_t currentStep )
        {
            typedef SuperCellDescription<
                SuperCellSize,
                typename CurlB::LowerMargin,
                typename CurlB::UpperMargin
            > BlockArea;
            AreaMapping< AREA, MappingDesc > mapper( cellDescription );
            constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                pmacc::math::CT::volume< SuperCellSize >::type::value
            >::value;
            // Note: here it is possible to check if PML is enabled in the local
            // domain at all, and otherwise optimize by calling the normal Yee update,
            // but we do not do that, as this is fragile wrt plans for close future changes
            PMACC_KERNEL( yeePML::KernelUpdateE< numWorkers, BlockArea >{ } )
                ( mapper.getGridDim(), numWorkers )(
                    CurlE(),
                    this->splitFields->getDeviceDataBox(),
                    this->fieldE->getDeviceDataBox(),
                    this->fieldB->getDeviceDataBox(),
                    mapper,
                    computeParameters( currentStep )
                );
        }

    public:

        using NummericalCellType = picongpu::numericalCellTypes::YeeCell;
        using CurrentInterpolation = T_CurrentInterpolation;
        using CurlE = T_CurlE;
        using CurlB = T_CurlB;

        YeePML( MappingDesc cellDescription ) :
            cellDescription( cellDescription )
        {
            // Split fields are created here to not waste memory in case
            // PML is not used
            DataConnector &dc = Environment<>::get().DataConnector();
            fieldE = dc.get< FieldE >( FieldE::getName(), true );
            fieldB = dc.get< FieldB >( FieldB::getName(), true );
            dc.share(std::shared_ptr< ISimulationData >(new yeePML::SplitFields(cellDescription)));
            splitFields = dc.get< yeePML::SplitFields >( yeePML::SplitFields::getName(), true );
            initializeParameters();
        }

        void update_beforeCurrent( uint32_t const currentStep )
        {
            // Steps are the same as in the Yee solver,
            // PML updates are done as part of updateE(), updateBHalf()
            updateBHalf < CORE + BORDER >( currentStep );
            EventTask eRfieldB = fieldB->asyncCommunication( __getTransactionEvent() );

            updateE< CORE >( currentStep );
            __setTransactionEvent( eRfieldB );
            updateE< BORDER >( currentStep );
        }

        void update_afterCurrent( uint32_t const currentStep )
        {
            // Steps are the same as in the Yee solver,
            // except the FieldManipulator::absorbBorder() is not called,
            // PML updates are done as part of updateE(), updateBHalf()
            if( laserProfiles::Selected::INIT_TIME > 0.0_X )
                LaserPhysics{}( currentStep );

            EventTask eRfieldE = fieldE->asyncCommunication( __getTransactionEvent() );

            updateBHalf< CORE >( currentStep );
            __setTransactionEvent( eRfieldE );
            updateBHalf< BORDER >( currentStep );

            EventTask eRfieldB = fieldB->asyncCommunication( __getTransactionEvent() );
            __setTransactionEvent( eRfieldB );
        }
    };

} // namespace maxwellSolver
} // namespace fields
} // namespace picongpu

#include "picongpu/fields/MaxwellSolver/YeePML/SplitFields.tpp"
