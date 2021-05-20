/* Copyright 2021 Sergei Bastrakov
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

#include "picongpu/fields/absorber/Absorber.hpp"

#include <boost/program_options/options_description.hpp>

#include <cstdint>
#include <stdexcept>
#include <string>


namespace picongpu
{
    namespace simulation
    {
        namespace stage
        {
            /** Functor for the stage of the PIC loop performing field absorption
             *
             * This stage does not run by itself, but is needed to propagate command-line parameters
             */
            class FieldAbsorber
            {
            public:
                /** Register program options for field absorber
                 *
                 * @param desc program options following boost::program_options::options_description
                 */
                void registerHelp(po::options_description& desc)
                {
                    desc.add_options()(
                        "fieldAbsorber",
                        po::value<std::string>(&kindName),
                        (std::string("Field absorber kind [exponential, pml] default: " + kindName).c_str()));
                }

                /** Initialize the field absorber stage
                 *
                 * This method has to be called during initialization of the simulation.
                 * Before this method is called, the instance of fields::absorber::Absorber cannot be used safely.
                 */
                void init()
                {
                    using namespace fields::absorber;
                    // So far there are only two kinds and so names are hardcoded
                    if(kindName == "exponential")
                        Absorber::kind() = Absorber::Kind::Exponential;
                    else if(kindName == "pml")
                        Absorber::kind() = Absorber::Kind::Pml;
                    else
                        throw std::runtime_error("Unsupported field absorber type");
                }

            private:
                //! Name set by program option
                std::string kindName = "exponential";
            };

        } // namespace stage
    } // namespace simulation
} // namespace picongpu
