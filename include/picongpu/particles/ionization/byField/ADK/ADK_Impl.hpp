/* Copyright 2015-2017 Marco Garten
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
#include <pmacc/traits/Resolve.hpp>
#include "picongpu/traits/UsesRNG.hpp"

#include "picongpu/fields/FieldB.hpp"
#include "picongpu/fields/FieldE.hpp"

#include "picongpu/particles/ionization/byField/ADK/ADK.def"
#include "picongpu/particles/ionization/byField/ADK/AlgorithmADK.hpp"

#include <pmacc/random/methods/AlpakaRand.hpp>
#include <pmacc/random/distributions/Uniform.hpp>
#include <pmacc/random/RNGProvider.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/compileTime/conversion/TypeToPointerPair.hpp>
#include <pmacc/memory/boxes/DataBox.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>

#include <boost/type_traits/integral_constant.hpp>


namespace picongpu
{
namespace traits
{
    /** specialization of the UsesRNG trait
     * --> ionization module uses random number generation
     */
    template<typename T_IonizationAlgorithm, typename T_DestSpecies, typename T_SrcSpecies>
    struct UsesRNG<particles::ionization::ADK_Impl<T_IonizationAlgorithm, T_DestSpecies, T_SrcSpecies> > :
    public boost::true_type
    {
    };
} // namespace traits

namespace particles
{
namespace ionization
{

    /** \struct ADK_Impl
     *
     * \brief Ammosov-Delone-Krainov
     *        Tunneling ionization for hydrogenlike atoms
     *
     * \tparam T_DestSpecies electron species to be created
     * \tparam T_SrcSpecies particle species that is ionized
     */
    template<typename T_IonizationAlgorithm, typename T_DestSpecies, typename T_SrcSpecies>
    struct ADK_Impl
    {

        typedef T_DestSpecies DestSpecies;
        typedef T_SrcSpecies  SrcSpecies;

        typedef typename SrcSpecies::FrameType FrameType;

        /* specify field to particle interpolation scheme */
        typedef typename pmacc::traits::Resolve<
            typename GetFlagType<FrameType,interpolation<> >::type
        >::type Field2ParticleInterpolation;

        /* margins around the supercell for the interpolation of the field on the cells */
        typedef typename GetMargin<Field2ParticleInterpolation>::LowerMargin LowerMargin;
        typedef typename GetMargin<Field2ParticleInterpolation>::UpperMargin UpperMargin;

        /* relevant area of a block */
        typedef SuperCellDescription<
            typename MappingDesc::SuperCellSize,
            LowerMargin,
            UpperMargin
            > BlockArea;

        BlockArea BlockDescription;

        private:

            /* define ionization ALGORITHM (calculation) for ionization MODEL */
            typedef T_IonizationAlgorithm IonizationAlgorithm;

            /* random number generator */
            typedef pmacc::random::RNGProvider<simDim, pmacc::random::methods::AlpakaRand< cupla::Acc>> RNGFactory;
            typedef pmacc::random::distributions::Uniform<float> Distribution;
            typedef typename RNGFactory::GetRandomType<Distribution>::type RandomGen;
            RandomGen randomGen;

            typedef MappingDesc::SuperCellSize TVec;

            typedef FieldE::ValueType ValueType_E;
            typedef FieldB::ValueType ValueType_B;
            /* global memory EM-field device databoxes */
            PMACC_ALIGN(eBox, FieldE::DataBoxType);
            PMACC_ALIGN(bBox, FieldB::DataBoxType);
            /* shared memory EM-field device databoxes */
            PMACC_ALIGN(cachedE, DataBox<SharedBox<ValueType_E, typename BlockArea::FullSuperCellSize,1> >);
            PMACC_ALIGN(cachedB, DataBox<SharedBox<ValueType_B, typename BlockArea::FullSuperCellSize,0> >);

        public:
            /* host constructor initializing member : random number generator */
            ADK_Impl(const uint32_t currentStep) : randomGen(RNGFactory::createRandom<Distribution>())
            {
                DataConnector &dc = Environment<>::get().DataConnector();
                /* initialize pointers on host-side E-(B-)field databoxes */
                auto fieldE = dc.get< FieldE >( FieldE::getName(), true );
                auto fieldB = dc.get< FieldB >( FieldB::getName(), true );
                /* initialize device-side E-(B-)field databoxes */
                eBox = fieldE->getDeviceDataBox();
                bBox = fieldB->getDeviceDataBox();

            }

            /** Initialization function on device
             *
             * \brief Cache EM-fields on device
             *         and initialize possible prerequisites for ionization, like e.g. random number generator.
             *
             * This function will be called inline on the device which must happen BEFORE threads diverge
             * during loop execution. The reason for this is the `__syncthreads()` call which is necessary after
             * initializing the E-/B-field shared boxes in shared memory.
             */
            template< typename T_Acc >
            DINLINE void init(
                T_Acc const & acc,
                const DataSpace<simDim>& blockCell,
                const int& linearThreadIdx,
                const DataSpace<simDim>& localCellOffset
            )
            {

                /* caching of E and B fields */
                cachedB = CachedBox::create<
                    0,
                    ValueType_B
                >(
                    acc,
                    BlockArea()
                );
                cachedE = CachedBox::create<
                    1,
                    ValueType_E
                >(
                    acc,
                    BlockArea()
                );

                /* instance of nvidia assignment operator */
                nvidia::functors::Assign assign;
                /* copy fields from global to shared */
                auto fieldBBlock = bBox.shift(blockCell);
                ThreadCollective<
                    BlockArea,
                    pmacc::math::CT::volume< typename BlockArea::SuperCellSize >::type::value
                > collective( linearThreadIdx );
                collective(
                          acc,
                          assign,
                          cachedB,
                          fieldBBlock
                          );
                /* copy fields from global to shared */
                auto fieldEBlock = eBox.shift(blockCell);
                collective(
                          acc,
                          assign,
                          cachedE,
                          fieldEBlock
                          );

                /* wait for shared memory to be initialized */
                __syncthreads();

                /* initialize random number generator with the local cell index in the simulation */
                this->randomGen.init(localCellOffset);
            }

            /** Determine number of new macro electrons due to ionization
             *
             * \param ionFrame reference to frame of the to-be-ionized particles
             * \param localIdx local (linear) index in super cell / frame
             */
            template< typename T_Acc >
            DINLINE uint32_t numNewParticles(const T_Acc& acc, FrameType& ionFrame, int localIdx)
            {
                /* alias for the single macro-particle */
                auto particle = ionFrame[localIdx];
                /* particle position, used for field-to-particle interpolation */
                floatD_X pos = particle[position_];
                const int particleCellIdx = particle[localCellIdx_];
                /* multi-dim coordinate of the local cell inside the super cell */
                DataSpace<TVec::dim> localCell(DataSpaceOperations<TVec::dim>::template map<TVec > (particleCellIdx));
                /* interpolation of E- */
                const fieldSolver::numericalCellType::traits::FieldPosition<FieldE> fieldPosE;
                ValueType_E eField = Field2ParticleInterpolation()
                    (cachedE.shift(localCell).toCursor(), pos, fieldPosE());
                /*                     and B-field on the particle position */
                const fieldSolver::numericalCellType::traits::FieldPosition<FieldB> fieldPosB;
                ValueType_B bField = Field2ParticleInterpolation()
                    (cachedB.shift(localCell).toCursor(), pos, fieldPosB());

                /* define number of bound macro electrons before ionization */
                float_X prevBoundElectrons = particle[boundElectrons_];

                IonizationAlgorithm ionizeAlgo;
                /* determine number of new macro electrons to be created */
                uint32_t newMacroElectrons = ionizeAlgo(
                                                bField, eField,
                                                particle, this->randomGen(acc)
                                                );

                return newMacroElectrons;

            }

            /* Functor implementation
             *
             * Ionization model specific particle creation
             *
             * \tparam T_parentIon type of ion species that is being ionized
             * \tparam T_childElectron type of electron species that is created
             * \param parentIon ion instance that is ionized
             * \param childElectron electron instance that is created
             */
            template<typename T_parentIon, typename T_childElectron, typename T_Acc>
            DINLINE void operator()(const T_Acc& acc, T_parentIon& parentIon,T_childElectron& childElectron)
            {
                /* for not mixing operations::assign up with the nvidia functor assign */
                namespace partOp = pmacc::particles::operations;
                /* each thread sets the multiMask hard on "particle" (=1) */
                childElectron[multiMask_] = 1u;
                const float_X weighting = parentIon[weighting_];

                /* each thread initializes a clone of the parent ion but leaving out
                 * some attributes:
                 * - multiMask: reading from global memory takes longer than just setting it again explicitly
                 * - momentum: because the electron would get a higher energy because of the ion mass
                 * - boundElectrons: because species other than ions or atoms do not have them
                 * (gets AUTOMATICALLY deselected because electrons do not have this attribute)
                 */
                auto targetElectronClone = partOp::deselect<bmpl::vector2<multiMask, momentum> >(childElectron);

                partOp::assign(targetElectronClone, partOp::deselect<particleId>(parentIon));

                const float_X massIon = attribute::getMass(weighting,parentIon);
                const float_X massElectron = attribute::getMass(weighting,childElectron);

                const float3_X electronMomentum (parentIon[momentum_]*(massElectron/massIon));

                childElectron[momentum_] = electronMomentum;

                /* conservation of momentum
                 * \todo add conservation of mass */
                parentIon[momentum_] -= electronMomentum;

                /** ionization of the ion by reducing the number of bound electrons
                 *
                 * @warning substracting a float from a float can potentially
                 *          create a negative boundElectrons number for the ion,
                 *          see #1850 for details
                 */
                parentIon[boundElectrons_] -= float_X(1.);
            }

    };

} // namespace ionization
} // namespace particles
} // namespace picongpu