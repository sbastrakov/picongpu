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
#include "picongpu/fields/MaxwellSolver/YeePML/Field.hpp"
#include "picongpu/fields/MaxwellSolver/YeePML/Parameters.hpp"
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

    /* Note: the yeePML namespace is only used for details and not the solver
     * itself in order to be consistent with other field solvers.
     */
    namespace yeePML
    {
    namespace detail
    {

        /**
         * Implementation of Yee + PML solver updates of E and B.
         *
         * @tparam T_CurlE functor to compute curl of E
         * @tparam T_CurlB functor to compute curl of B
         */
        template<
            typename T_CurlE,
            typename T_CurlB
        >
        class Solver
        {
        public:

            using CurlE = T_CurlE;
            using CurlB = T_CurlB;

            // Helper types for configuring kernels
            template< typename T_Curl >
            using BlockDescription = pmacc::SuperCellDescription<
                SuperCellSize,
                typename T_Curl::LowerMargin,
                typename T_Curl::UpperMargin
            >;
            template< uint32_t T_Area >
            using AreaMapper = pmacc::AreaMapping<
                T_Area,
                MappingDesc
            >;

            Solver( MappingDesc cellDescription ) :
                cellDescription{ cellDescription }
            {
                initParameters( );
                initFields( );
            }

            //! Get a reference to (full) field E
            picongpu::FieldE& getFieldE( )
            {
                return *( fieldE.get() );
            }

            //! Get a reference to (full) field B
            picongpu::FieldB& getFieldB( )
            {
                return *( fieldB.get() );
            }

            /**
             * Propagate B values in the given area by half a time step.
             *
             * @tparam T_Area area to apply updates to, the curl must be
             * applicable to all points; normally CORE, BORDER, or CORE + BORDER
             *
             * @param currentStep index of the current time iteration
             */
            template< uint32_t T_Area >
            void updateBHalf( uint32_t const currentStep )
            {
                constexpr auto numWorkers = getNumWorkers();
                using Kernel = yeePML::KernelUpdateBHalf<
                    numWorkers,
                    BlockDescription< CurlE >
                >;
                AreaMapper< T_Area > mapper{ cellDescription };
                /* Note: here it is possible to first check if PML is enabled
                 * in the local domain at all, and otherwise optimize by calling
                 * the normal Yee update kernel. We do not do that, as this
                 * would be fragile wrt future separation of PML into a plugin.
                 */
                PMACC_KERNEL( Kernel{ } )
                    ( mapper.getGridDim(), numWorkers )(
                        CurlE( ),
                        splitB->getDeviceDataBox(),
                        fieldB->getDeviceDataBox(),
                        fieldE->getDeviceDataBox(),
                        mapper,
                        getLocalParameters( currentStep )
                    );
            }

            /**
            * Propagate E values in the given area by a time step.
            *
            * @tparam T_Area area to apply updates to, the curl must be
            * applicable to all points; normally CORE, BORDER, or CORE + BORDER
            *
            * @param currentStep index of the current time iteration
            */
            template< uint32_t T_Area >
            void updateE( uint32_t currentStep )
            {
                constexpr auto numWorkers = getNumWorkers();
                using Kernel = yeePML::KernelUpdateE<
                    numWorkers,
                    BlockDescription< CurlB >
                >;
                AreaMapper< T_Area > mapper{ cellDescription };
                // Note: optimization considerations same as in updateBHalf().
                PMACC_KERNEL( Kernel{ } )
                    ( mapper.getGridDim(), numWorkers )(
                        CurlB(),
                        splitE->getDeviceDataBox(),
                        fieldE->getDeviceDataBox(),
                        fieldB->getDeviceDataBox(),
                        mapper,
                        getLocalParameters( currentStep )
                    );
            }

        private:

            // Yee solver data
            std::shared_ptr< picongpu::FieldE > fieldE;
            std::shared_ptr< picongpu::FieldB > fieldB;
            MappingDesc cellDescription;

            // PML field data
            std::shared_ptr< yeePML::FieldE > splitE;
            std::shared_ptr< yeePML::FieldB > splitB;

            // PML parameters
            /**
            * Thickness in terms of the global domain.
            *
            * We store only global thickness, as the local one can change
            * during the simulation and so has to be recomputed each time step.
            * There are no limitations on the size, as long as it fits a single
            * layer of external local domains (near the global simulation area
            * boundary) on each time step. In particular, PML is not required
            * to be aligned with the BORDER area.
            */
            Thickness globalSize;
            Parameters parameters;

            void initParameters( )
            {
                globalSize = getGlobalThickness( );
                parameters.sigmaKappaGradingOrder = SIGMA_GRADING_ORDER;
                parameters.alphaGradingOrder = ALPHA_GRADING_ORDER;
                for( auto dim = 0; dim < simDim; dim++ )
                {
                    parameters.normalizedSigmaMax[ dim ] = NORMALIZED_SIGMA_MAX[ dim ];
                    parameters.kappaMax[ dim ] = KAPPA_MAX[ dim ];
                    parameters.normalizedAlphaMax[ dim ] = NORMALIZED_ALPHA_MAX[ dim ];
                }
            }

            Thickness getGlobalThickness( ) const
            {
                Thickness globalThickness;
                for( auto axis = 0; axis < simDim; axis++ )
                    for( auto direction = 0; direction < 2; direction++ )
                        globalThickness( axis, direction ) = NUM_CELLS[ axis ][ direction ];
                return globalThickness;
            }

            void initFields()
            {
                /* Split fields are created here (and not with normal E and B)
                * in order to not waste memory in case PML is not used.
                */
                DataConnector & dc = Environment<>::get().DataConnector();
                fieldE = dc.get< picongpu::FieldE >( picongpu::FieldE::getName(), true );
                fieldB = dc.get< picongpu::FieldB >( picongpu::FieldB::getName(), true );
                dc.share(
                    std::shared_ptr< ISimulationData >(
                        new yeePML::FieldE{ cellDescription }
                        )
                );
                dc.share(
                    std::shared_ptr< ISimulationData >(
                        new yeePML::FieldB{ cellDescription }
                        )
                );
                splitE = dc.get< yeePML::FieldE >( yeePML::FieldE::getName(), true );
                splitB = dc.get< yeePML::FieldB >( yeePML::FieldB::getName(), true );
            }

            yeePML::detail::LocalParameters getLocalParameters(
                uint32_t const currentStep
            ) const
            {
                Thickness localThickness = getLocalThickness( currentStep );
                checkLocalThickness( localThickness );
                return yeePML::detail::LocalParameters( parameters, localThickness );
            }

            /**
             * Get PML thickness for the local domain at the current time step.
             * It depends on the current step because of the moving window.
             */
            Thickness getLocalThickness( uint32_t const currentStep ) const
            {
                /* The logic of the following checks is the same as in
                 * FieldManipulator::absorbBorder(), to disable the absorber
                 * at a border we set the corresponding thickness to 0.
                 */
                auto & movingWindow = MovingWindow::getInstance( );
                auto const numSlides = movingWindow.getSlideCounter( currentStep );
                auto const numExchanges = NumberOfExchanges< simDim >::value;
                auto const communicationMask = Environment< simDim >::get( ).GridController( ).getCommunicationMask( );
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

            //! Verify that PML fits the local domain
            void checkLocalThickness( Thickness const localThickness ) const
            {
                auto const localDomain = Environment< simDim >::get( ).SubGrid( ).getLocalDomain( );
                auto const localPMLSize = localThickness.negativeBorder + localThickness.positiveBorder;
                auto pmlFitsDomain = true;
                for( auto dim = 0; dim < simDim; dim++ )
                    if( localPMLSize[ dim ] > localDomain.size[ dim ] )
                        pmlFitsDomain = false;
                if( !pmlFitsDomain )
                    throw std::out_of_range("Requisted PML size exceeds the local domain");
            }

            //! Get number of workers for kernels
            static constexpr uint32_t getNumWorkers()
            {
                return pmacc::traits::GetNumWorkers<
                    pmacc::math::CT::volume< SuperCellSize >::type::value
                >::value;
            }

        };

    } // namespace detail
    } // namespace yeePML

    /**
      * Yee solver with perfectly matched layer (PML) absorber.
      *
      * Follows implicitly defined interface of field solvers.
      * Most of the actual solver is implemented by yeePML::detail::Solver.
      *
      * @tparam T_CurrentInterpolation current interpolation functor
      * (not used in the implementation, but part of field solver interface)
      * @tparam T_CurlE functor to compute curl of E
      * @tparam T_CurlB functor to compute curl of B
      */
    template<
        typename T_CurrentInterpolation,
        typename T_CurlE,
        typename T_CurlB
    >
    class YeePML
    {
    public:

        // Types required by field solver interface
        using NummericalCellType = picongpu::numericalCellTypes::YeeCell;
        using CurrentInterpolation = T_CurrentInterpolation;
        using CurlE = T_CurlE;
        using CurlB = T_CurlB;

        YeePML( MappingDesc cellDescription ) :
            solver( cellDescription )
        {
        }

        /**
         * Perform the first part of E and B propagation by a time step.
         *
         * Together with update_afterCurrent() forms the full propagation.
         *
         * @param currentStep index of the current time iteration
         */
        void update_beforeCurrent( uint32_t const currentStep )
        {
            /* These steps are the same as in the Yee solver,
             * PML updates are done as part of solver.updateE(),
             * solver.updateBHalf()
             */
            solver.updateBHalf < CORE + BORDER >( currentStep );
            auto & fieldB = solver.getFieldB( );
            EventTask eRfieldB = fieldB.asyncCommunication( __getTransactionEvent( ) );

            solver.updateE< CORE >( currentStep );
            __setTransactionEvent( eRfieldB );
            solver.updateE< BORDER >( currentStep );
        }

        /**
        * Perform the last part of E and B propagation by a time step.
        *
        * Together with update_beforeCurrent() forms the full propagation.
        *
        * @param currentStep index of the current time iteration
        */
        void update_afterCurrent( uint32_t const currentStep )
        {
            /* These steps are the same as in the Yee solver,
             * except the FieldManipulator::absorbBorder() is not called,
             * PML updates are done as part of solver.updateBHalf().
             */
            if( laserProfiles::Selected::INIT_TIME > 0.0_X )
                LaserPhysics{ }( currentStep );

            auto & fieldE = solver.getFieldE( );
            EventTask eRfieldE = fieldE.asyncCommunication( __getTransactionEvent( ) );

            solver.updateBHalf< CORE >( currentStep );
            __setTransactionEvent( eRfieldE );
            solver.updateBHalf< BORDER >( currentStep );

            auto & fieldB = solver.getFieldB();
            EventTask eRfieldB = fieldB.asyncCommunication( __getTransactionEvent( ) );
            __setTransactionEvent( eRfieldB );
        }

        static pmacc::traits::StringProperty getStringProperties( )
        {
            pmacc::traits::StringProperty propList( "name", "YeePML" );
            return propList;
        }

    private:

        yeePML::detail::Solver< CurlE, CurlB > solver;

    };

} // namespace maxwellSolver
} // namespace fields
} // namespace picongpu

#include "picongpu/fields/MaxwellSolver/YeePML/Field.tpp"
