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

#include "picongpu/simulation/stage/FieldBackground.hpp"


namespace picongpu
{
    namespace simulation
    {
        namespace stage
        {


                void FieldBackground::registerHelp(po::options_description& desc)
                {
                    desc.add_options()(
                        "fieldBackground.duplicateFields",
                        po::value<bool>(&duplicateFields)->zero_tokens(),
                        "duplicate E and B field storage inside field background to improve its performance "
                        "and potentially avoid some numerical noise at cost of using more memory, "
                        "only affects the fields with activated background");
                }

                /** Initialize field background stage
                 *
                 * This method must be called once before calling add(), subtract() and fillSimulation().
                 * The initialization has to be delayed for this class as it needs registerHelp() like the plugins do.
                 *
                 * @param cellDescription mapping for kernels
                 */
                void FieldBackground::init(MappingDesc const cellDescription)
                {
                    applyE = std::make_unique<ApplyE>(cellDescription, duplicateFields);
                    applyB = std::make_unique<ApplyB>(cellDescription, duplicateFields);
                }

                /** Add field background to the electromagnetic field
                 *
                 * Affects data sets named FieldE::getName(), FieldB::getName().
                 * As the result of this operation, they will have a sum of old values and background values.
                 *
                 * @param step index of time iteration
                 */
                void FieldBackground::add(uint32_t const step)
                {
                    checkInitialization();
                    applyE->add(step);
                    applyB->add(step);
                }

                /** Subtract field background from the electromagnetic field
                 *
                 * Affects data sets named FieldE::getName(), FieldB::getName().
                 * As the result of this operation, they will have values like before the last call to add().
                 *
                 * Warning: when fieldBackground.duplicateFields is enabled, the fields are assumed to not have changed
                 * since the call to add(). Having fieldBackground.duplicateFields disabled does not rely on this.
                 * However, this assumption should generally hold true in the PIC computational loop.
                 *
                 * @param step index of time iteration
                 */
                void FieldBackground::subtract(uint32_t const step)
                {
                    checkInitialization();
                    applyE->subtract(step);
                    applyB->subtract(step);
                }

                //! Check if this class was properly initialized, throws when failed
                void FieldBackground::checkInitialization() const
                {
                    if(!applyE || !applyB)
                        throw std::runtime_error("simulation::stage::FieldBackground used without init() called");
                }

        } // namespace stage
    } // namespace simulation
} // namespace picongpu
