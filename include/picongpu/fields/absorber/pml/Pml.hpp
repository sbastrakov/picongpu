/* Copyright 2019-2021 Sergei Bastrakov
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

#include "picongpu/fields/absorber/Absorber.hpp"
#include "picongpu/fields/absorber/pml/Field.hpp"
#include "picongpu/fields/absorber/pml/Parameters.hpp"
#include "picongpu/fields/absorber/pml/Pml.kernel"

#include <cstdint>
#include <string>


namespace picongpu
{
    namespace fields
    {
        namespace absorber
        {
            namespace pml
            {
                /** Implementation of Yee + PML solver updates of E and B
                 *
                 * The original paper on this approach is J.A. Roden, S.D. Gedney.
                 * Convolution PML (CPML): An efficient FDTD implementation of the
                 * CFS - PML for arbitrary media. Microwave and optical technology
                 * letters. 27 (5), 334-339 (2000).
                 * https://doi.org/10.1002/1098-2760(20001205)27:5%3C334::AID-MOP14%3E3.0.CO;2-A
                 * Our implementation is based on a more detailed description in section
                 * 7.9 of the book A. Taflove, S.C. Hagness. Computational
                 * Electrodynamics. The Finite-Difference Time-Domain Method. Third
                 * Edition. Artech house, Boston (2005), referred to as
                 * [Taflove, Hagness].
                 */
                class Pml : public Absorber
                {
                public:
                    //! Create PML absorber instance
                    Pml()
                    {
                        // Copy thickness from pml.param
                        for(uint32_t axis = 0u; axis < 3u; axis++)
                            for(uint32_t direction = 0u; direction < 2u; direction++)
                                numCells[axis][direction] = maxwellSolver::Pml::NUM_CELLS[axis][direction];
                        name = std::string{"convolutional PML"};
                        initParameters();
                    }

                    /** Initialize PML field absorber
                     *
                     * @param newCellDescription mapping for kernels
                     */
                    void init(MappingDesc const newCellDescription) override
                    {
                        Absorber::init(newCellDescription);
                        DataConnector& dc = Environment<>::get().DataConnector();
                        psiE = std::make_shared<pml::FieldE>(cellDescription, getGlobalThickness());
                        psiB = std::make_shared<pml::FieldB>(cellDescription, getGlobalThickness());
                        dc.share(psiE);
                        dc.share(psiB);
                    }

                    /** Functor to update electric field by a time step using FDTD with the given curl and PML
                     *
                     * @tparam T_CurlB curl functor type according to the Curl concept
                     * @tparam T_Area area to apply updates to
                     *
                     * @param currentStep index of the current time iteration
                     */
                    template<typename T_CurlB, uint32_t T_Area>
                    UpdateEFunctor<T_CurlB> getUpdateEFunctor(uint32_t const currentStep)
                    {
                        AreaMapper<T_Area> mapper{cellDescription};
                        return UpdateEFunctor<T_CurlB>{
                            psiE->getDeviceOuterLayerBox(),
                            getLocalParameters(mapper, currentStep)};
                    }

                    /** Functor to update magnetic field by half a time step using FDTD with the given curl and PML
                     *
                     * @tparam T_CurlE curl functor type according to the Curl concept
                     * @tparam T_Area area to apply updates to
                     *
                     * @param currentStep index of the current time iteration
                     * @param updatePsiB whether convolutional magnetic fields need to be updated, or are
                     * up-to-date
                     */
                    template<typename T_CurlE, uint32_t T_Area>
                    UpdateBHalfFunctor<T_CurlE> getUpdateBHalfFunctor(
                        uint32_t const currentStep,
                        bool const updatePsiB)
                    {
                        AreaMapper<T_Area> mapper{cellDescription};
                        return UpdateBHalfFunctor<T_CurlE>{
                            psiB->getDeviceOuterLayerBox(),
                            getLocalParameters(mapper, currentStep),
                            updatePsiB};
                    }

                private:
                    void initParameters()
                    {
                        namespace paramPml = maxwellSolver::Pml;
                        parameters.sigmaKappaGradingOrder = paramPml::SIGMA_KAPPA_GRADING_ORDER;
                        parameters.alphaGradingOrder = paramPml::ALPHA_GRADING_ORDER;
                        for(uint32_t dim = 0u; dim < simDim; dim++)
                        {
                            parameters.normalizedSigmaMax[dim] = paramPml::NORMALIZED_SIGMA_MAX[dim];
                            parameters.kappaMax[dim] = paramPml::KAPPA_MAX[dim];
                            parameters.normalizedAlphaMax[dim] = paramPml::NORMALIZED_ALPHA_MAX[dim];
                        }
                    }


                    template<uint32_t T_Area>
                    using AreaMapper = pmacc::AreaMapping<T_Area, MappingDesc>;

                    template<uint32_t T_Area>
                    LocalParameters getLocalParameters(AreaMapper<T_Area>& mapper, uint32_t const currentStep) const
                    {
                        Thickness localThickness = getLocalThickness(currentStep);
                        checkLocalThickness(localThickness);
                        return LocalParameters(
                            parameters,
                            localThickness,
                            mapper.getGridSuperCells() * SuperCellSize::toRT(),
                            mapper.getGuardingSuperCells() * SuperCellSize::toRT());
                    }

                    /** Get PML thickness for the local domain at the current time step.
                     *
                     * It depends on the current step because of the moving window.
                     */
                    Thickness getLocalThickness(uint32_t const currentStep) const
                    {
                        /* The logic of the following checks is to disable the absorber
                         * at a border we set the corresponding thickness to 0.
                         */
                        auto& movingWindow = MovingWindow::getInstance();
                        auto const numSlides = movingWindow.getSlideCounter(currentStep);
                        auto const numExchanges = NumberOfExchanges<simDim>::value;
                        auto const communicationMask
                            = Environment<simDim>::get().GridController().getCommunicationMask();
                        Thickness localThickness = getGlobalThickness();
                        for(uint32_t exchange = 1u; exchange < numExchanges; ++exchange)
                        {
                            /* Here we are only interested in the positive and negative
                             * directions for x, y, z axes and not the "diagonal" ones.
                             * So skip other directions except left, right, top, bottom,
                             * back, front
                             */
                            if(FRONT % exchange != 0)
                                continue;

                            // Transform exchange into a pair of axis and direction
                            uint32_t axis = 0;
                            if(exchange >= BOTTOM && exchange <= TOP)
                                axis = 1;
                            if(exchange >= BACK)
                                axis = 2;
                            uint32_t direction = exchange % 2;

                            // No PML at the borders between two local domains
                            bool hasNeighbour = communicationMask.isSet(exchange);
                            if(hasNeighbour)
                                localThickness(axis, direction) = 0;

                            // Disable PML during laser initialization
                            if(fields::laserProfiles::Selected::initPlaneY == 0)
                            {
                                bool isLaserInitializationOver
                                    = (currentStep * DELTA_T) >= fields::laserProfiles::Selected::INIT_TIME;
                                if(numSlides == 0 && !isLaserInitializationOver && exchange == TOP)
                                    localThickness(axis, direction) = 0;
                            }

                            // Disable PML at the far side of the moving window
                            if(movingWindow.isSlidingWindowActive(currentStep) && exchange == BOTTOM)
                                localThickness(axis, direction) = 0;
                        }
                        return localThickness;
                    }

                    //! Verify that PML fits the local domain
                    void checkLocalThickness(Thickness const localThickness) const
                    {
                        auto const localDomain = Environment<simDim>::get().SubGrid().getLocalDomain();
                        auto const localPMLSize
                            = localThickness.getNegativeBorder() + localThickness.getPositiveBorder();
                        auto pmlFitsDomain = true;
                        for(uint32_t dim = 0u; dim < simDim; dim++)
                            if(localPMLSize[dim] > localDomain.size[dim])
                                pmlFitsDomain = false;
                        if(!pmlFitsDomain)
                            throw std::out_of_range("Requested PML size exceeds the local domain");
                    }

                    /// TODO update the comment
                    /** Thickness in terms of the global domain.
                     *
                     * We store only global thickness, as the local one can change
                     * during the simulation and so has to be recomputed for each time
                     * step. PML must be fully contained in a single layer of local
                     * domains near the global simulation area boundary. (Note that
                     * the domains of this layer might be changing, e.g. due to moving
                     * window.) There are no other limitations on PML thickness. In
                     * particular, it is independent of the BORDER area size.
                     */
                    Parameters parameters;

                    /* PML convolutional field data, defined as in [Taflove, Hagness],
                     * eq. (7.105a,b), and similar for other components
                     */
                    std::shared_ptr<pml::FieldE> psiE;
                    std::shared_ptr<pml::FieldB> psiB;
                };

            } // namespace pml
        } // namespace absorber
    } // namespace fields
} // namespace picongpu

#include "picongpu/fields/absorber/pml/Field.tpp"
