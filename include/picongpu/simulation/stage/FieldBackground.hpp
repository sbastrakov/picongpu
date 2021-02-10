/* Copyright 2013-2021 Axel Huebl, Felix Schmitt, Heiko Burau, Rene Widera,
 *                     Richard Pausch, Alexander Debus, Marco Garten,
 *                     Benjamin Worpitz, Alexander Grund, Sergei Bastrakov
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

#include "picongpu/fields/background/cellwiseOperation.hpp"
#include "picongpu/fields/FieldB.hpp"
#include "picongpu/fields/FieldE.hpp"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/Environment.hpp>
#include <pmacc/nvidia/functors/Add.hpp>
#include <pmacc/nvidia/functors/Sub.hpp>
#include <pmacc/type/Area.hpp>

#include <cstdint>


namespace picongpu
{
    namespace simulation
    {
        namespace stage
        {
            //! Functor for the stage of the PIC loop applying field background
            class FieldBackground
            {
            public:
                /** Create field background functor with a fake mapping description
                 *
                 * The actual field background has to be set later by calling setMappingDescription()
                 */
                FieldBackground() : cellDescription(DataSpace<simDim>(SuperCellSize::toRT()))
                {
                }

                /** Set mapping description for kernels
                 *
                 * @param newCellDescription new mapping
                 */
                void setMappingDescription(MappingDesc const newCellDescription)
                {
                    cellDescription = newCellDescription;
                }

                /** Apply the field background to the electromagnetic field
                 *
                 * @param step index of time iteration
                 */
                void apply(uint32_t const step) const
                {
                    process<CORE + BORDER + GUARD>(step, pmacc::nvidia::functors::Add{});
                }

                /** Restore the original electromagnetic field by reverting application of the field background
                 *
                 * @param step index of time iteration
                 */
                void restore(uint32_t const step) const
                {
                    process<CORE + BORDER + GUARD>(step, pmacc::nvidia::functors::Sub{});
                }

                /** Apply initial field background during filling the simulation
                 *
                 * @param step index of time iteration
                 */
                void fillSimulation(uint32_t const step) const
                {
                    /* restore background fields in GUARD
                     *
                     * loads the outer GUARDS of the global domain for absorbing/open boundary condtions
                     *
                     * @todo as soon as we add GUARD fields to the checkpoint data, e.g. for PML boundary
                     *       conditions, this section needs to be removed
                     */
                    process<GUARD>(step, pmacc::nvidia::functors::Add());
                }

            private:
                /** Apply the given functor to the field background in the given area
                 *
                 * @tparam T_area area to operate on
                 * @tparam T_Functor functor type compatible to pmacc::nvidia::functors
                 *
                 * @param step index of time iteration
                 * @param functor functor to apply
                 */
                template<uint32_t T_area, typename T_Functor>
                void process(uint32_t const step, T_Functor functor) const
                {
                    using namespace pmacc;
                    DataConnector& dc = Environment<>::get().DataConnector();
                    auto fieldE = dc.get<FieldE>(FieldE::getName(), true);
                    auto fieldB = dc.get<FieldB>(FieldB::getName(), true);
                    using Background = cellwiseOperation::CellwiseOperation<T_area>;
                    Background background(cellDescription);
                    background(
                        fieldE,
                        functor,
                        FieldBackgroundE(fieldE->getUnit()),
                        step,
                        FieldBackgroundE::InfluenceParticlePusher);
                    background(
                        fieldB,
                        functor,
                        FieldBackgroundB(fieldB->getUnit()),
                        step,
                        FieldBackgroundB::InfluenceParticlePusher);
                    dc.releaseData(FieldE::getName());
                    dc.releaseData(FieldB::getName());
                }

                //! Mapping for kernels
                MappingDesc cellDescription;
            };

        } // namespace stage
    } // namespace simulation
} // namespace picongpu
