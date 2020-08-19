/* Copyright 2013-2020 Axel Huebl, Rene Widera
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

#include "picongpu/traits/SplashToPIC.hpp"
#include "picongpu/traits/PICToSplash.hpp"

#include "picongpu/plugins/PhaseSpace/AxisDescription.hpp"
#include <pmacc/communication/manager_common.hpp>
#include <pmacc/mappings/simulation/GridController.hpp>
#include <pmacc/mappings/simulation/SubGrid.hpp>
#include <pmacc/dimensions/DataSpace.hpp>
#include <pmacc/cuSTL/container/HostBuffer.hpp>
#include <pmacc/math/vector/Int.hpp>
#include <pmacc/verify.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <utility>
#include <mpi.h>
///#include <splash/splash.h>

namespace picongpu
{
    class DumpHBuffer
    {
    private:
       typedef typename MappingDesc::SuperCellSize SuperCellSize;

    public:
        /** Dump the PhaseSpace host Buffer
         *
         * \tparam Type the HBuffers element type
         * \tparam int the HBuffers dimension
         * \param hBuffer const reference to the hBuffer, including guard cells in spatial dimension
         * \param axis_element plot to create: e.g. py, x from momentum/spatial-coordinate
         * \param unit sim unit of the buffer
         * \param strSpecies unique short hand name of the species
         * \param currentStep current time step
         * \param mpiComm communicator of the participating ranks
         */
        template<typename T_Type, int T_bufDim>
        void operator()( const pmacc::container::HostBuffer<T_Type, T_bufDim>& hBuffer,
                         const AxisDescription axis_element,
                         const std::pair<float_X, float_X> axis_p_range,
                         const float_64 pRange_unit,
                         const float_64 unit,
                         const std::string strSpecies,
                         const uint32_t currentStep,
                         MPI_Comm mpiComm ) const
        {

        }
    };

} /* namespace picongpu */
