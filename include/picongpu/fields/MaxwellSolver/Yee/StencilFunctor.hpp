/* Copyright 2021 Rene Widera, Sergei Bastrakov
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

#include <cstdint>


namespace picongpu
{
    namespace fields
    {
        namespace maxwellSolver
        {
            namespace yee
            {                
                /** Base stencil functor to update fields inside the kernel
                 *
                 * This class serves to define the interface requirements for stencil functor implementations.
                 * So if roughly defines a "concept".
                 *
                 * @tparam T_CurlB curl functor type to be applied to magnetic field,
                 *                 adheres to the Curl concept
                 */
                template<typename T_Curl>
                struct StencilFunctor
                {
                public:
                    //! Stencil requirements for lower margins
                    using LowerMargin = typename traits::GetLowerMargin<T_CurlB>::type;
                    //! Stencil requirements for upper margins
                    using UpperMargin = typename traits::GetUpperMargin<T_CurlB>::type;

                    /** Update field at the given position
                     *
                     * @tparam T_LocalBBox local magnetic field box type
                     * @tparam T_LocalEBox local electric field box type
                     *
                     * @param gridIndex index of the updated field element, with guards
                     * @param localB magnetic field box shifted to position gridIndex,
                     *               note that it is the box, not the value
                     * @param localE electric field box shifted to position gridIndex,
                     *               note that it is the box, not the value
                     *
                     * @return update the value pointed to by localE
                     */
                    template<typename T_LocalBBox, typename T_LocalEBox>
                    DINLINE void operator()(
                        pmacc::DataSpace<simDim> const& gridIndex,
                        T_LocalBBox const localB,
                        T_LocalEBox localE)
                    {
                        // Not implemented in this class, the child classes should implement it
                    }
                };

            } // namespace yee
        } // namespace maxwellSolver
    } // namespace fields
} // namespace picongpu
