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

#include "picongpu/simulation_defines.hpp"
#include "picongpu/fields/background/cellwiseOperation.hpp"
#include "picongpu/fields/FieldB.hpp"
#include "picongpu/fields/FieldE.hpp"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/Environment.hpp>
#include <pmacc/nvidia/functors/Add.hpp>
#include <pmacc/nvidia/functors/Sub.hpp>
#include <pmacc/type/Area.hpp>

#include <boost/program_options.hpp>

#include <cstdint>
#include <memory>
#include <stdexcept>


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
                /** Register program options for field background
                 *
                 * @param desc program options following boost::program_options::options_description
                 */
                void registerHelp(po::options_description& desc);

                /** Initialize field background stage
                 *
                 * This method must be called once before calling add(), subtract() and fillSimulation().
                 * The initialization has to be delayed for this class as it needs registerHelp() like the plugins do.
                 *
                 * @param cellDescription mapping for kernels
                 */
                void init(MappingDesc const cellDescription);

                /** Add field background to the electromagnetic field
                 *
                 * Affects data sets named FieldE::getName(), FieldB::getName().
                 * As the result of this operation, they will have a sum of old values and background values.
                 *
                 * @param step index of time iteration
                 */
                void add(uint32_t const step);

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
                void subtract(uint32_t const step);

            private:
                //! Check if this class was properly initialized, throws when failed
                void checkInitialization() const;

                //! Implememtation type to apply background to field E
                using ApplyE = detail::ApplyFieldBackground<FieldE, FieldBackgroundE>;

                //! Object to apply background to field E
                std::unique_ptr<ApplyE> applyE;

                //! Implememtation type to apply background to field B
                using ApplyB = detail::ApplyFieldBackground<FieldB, FieldBackgroundB>;

                //! Object to apply background to field B
                std::unique_ptr<ApplyB> applyB;

                //! Flag to store duplicates fields with enabled backgrounds
                bool duplicateFields = false;
            };

        } // namespace stage
    } // namespace simulation
} // namespace picongpu
