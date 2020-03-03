/* Copyright 2013-2020 Axel Huebl, Heiko Burau, Rene Widera, Richard Pausch,
 *                     Klaus Steiniger, Felix Schmitt, Benjamin Worpitz,
 *                     Juncheng E, Sergei Bastrakov
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

#include "picongpu/particles/traits/SpeciesEligibleForSolver.hpp"
#include "picongpu/plugins/ISimulationPlugin.hpp"
#include "picongpu/plugins/common/stringHelpers.hpp"
#include "picongpu/plugins/xrayDiffraction/ComputeDiffraction.hpp"
#include "picongpu/plugins/xrayDiffraction/ReciprocalSpace.hpp"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/dimensions/DataSpaceOperations.hpp>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>


namespace picongpu
{
namespace plugins
{
namespace xrayDiffraction
{

    using namespace pmacc;

    namespace po = boost::program_options;

    /** X-ray diffraction plugin
     * This SAXS plugin simulates the SAXS scattering pattern on-the-fly
     * from the particle positions obtained from PIConGPU.
     **/
    template< typename T_Species >
    class XrayDiffraction : public ISimulationPlugin
    {
    private:

        using ComputeDiffraction = detail::ComputeDiffraction< T_Species >;
        std::unique_ptr< ComputeDiffraction > pImpl;

        MappingDesc cellDescription;
        std::string notifyPeriod; 
        std::string prefix;

        /** Range of scattering vector
         * The scattering vector here is defined as
         * 4*pi*sin(theta)/lambda, where 2theta is the angle
         * between scattered and incident beam
         **/
        float3_X q_min, q_max;

        //! Number of scattering vectors
        DataSpace< 3 > numVectors;

    public:

        XrayDiffraction():
            prefix( T_Species::FrameType::getName() + std::string("_xrayDiffraction") ),
            cellDescription( DataSpace< simDim >::create( 0 ) )
        {
            Environment<>::get().PluginConnector().registerPlugin( this );
        }

        /**
         * This function represents what is actually calculated if the plugin
         * is called. Here, one only sets the particles pointer to the data of
         * the latest time step and calls the 'calculateSAXS'
         * function if for the actual time step radiation is to be calculated.
         * @param currentStep
         */
        void notify(uint32_t currentStep)
        {
            pImpl->operator()(
                currentStep,
                cellDescription
            );
        }

        void pluginRegisterHelp( po::options_description &desc ) override
        {
            desc.add_options()((prefix + ".period").c_str(),
                po::value<std::string>(&notifyPeriod),
                "enable plugin [for each n-th step]")(
                (prefix + ".qx_max").c_str(),
                po::value<float_X>(&q_max[0])->default_value(5),
                "reciprocal space range qx_max (A^-1)")(
                (prefix + ".qy_max").c_str(),
                po::value<float_X>(&q_max[1])->default_value(5),
                "reciprocal space range qy_max (A^-1)")(
                (prefix + ".qz_max").c_str(),
                po::value<float_X>(&q_max[2])->default_value(5),
                "reciprocal space range qz_max (A^-1)")(
                (prefix + ".qx_min").c_str(),
                po::value<float_X>(&q_min[0])->default_value(-5),
                "reciprocal space range qx_min (A^-1)")(
                (prefix + ".qy_min").c_str(),
                po::value<float_X>(&q_min[1])->default_value(-5),
                "reciprocal space range qy_min (A^-1)")(
                (prefix + ".qz_min").c_str(),
                po::value<float_X>(&q_min[2])->default_value(-5),
                "reciprocal space range qz_min (A^-1)")(
                (prefix + ".n_qx").c_str(),
                po::value<int>(&numVectors[0])->default_value(100),
                "Number of qx")((prefix + ".n_qy").c_str(),
                po::value<int>(&numVectors[1])->default_value(100),
                "Number of qy")((prefix + ".n_qz").c_str(),
                po::value<int>(&numVectors[2])->default_value(1), "Number of qz");
        }

        std::string pluginGetName() const override
        {
            return std::string{
                "X-ray diffraction: calculate diffraction scattering "
                "intensity of a species"
            };
        }

        void setMappingDescription( MappingDesc * newCellDescription ) override
        {
            cellDescription = *newCellDescription;
        }

        void restart(
            uint32_t restartStep,
            std::string const restartDirectory
        ) override
        {
            // No state to be read
        }

        void checkpoint(
            uint32_t currentStep,
            std::string const restartDirectory
        ) override
        {
            // No state to be saved
        }

    private:

        /**
         * The plugin is loaded on every MPI rank, and therefore this function is
         * executed on every MPI rank.
         * One host with MPI rank 0 is defined to be the master.
         * It creates a folder where all the
         * results are saved in a plain text format.
         **/
        void pluginLoad()
        {
            if (!notifyPeriod.empty())
            {
                auto qStep = ( q_max - q_min ) / precisionCast< float_X >( numVectors );
                auto reciprocalSpace = detail::ReciprocalSpace{
                    q_min,
                    qStep,
                    numVectors
                };
                pImpl = memory::makeUnique< ComputeDiffraction >(
                    reciprocalSpace,
                    prefix
                );
                Environment<>::get().PluginConnector().setNotificationPeriod(
                    this,
                    notifyPeriod
                );
            }
        }

        void pluginUnload() override
        {
            if( !notifyPeriod.empty() )
                CUDA_CHECK( cudaGetLastError() );
        }

    };

} // namespace xrayDiffraction
} // namespace plugins

namespace particles
{
namespace traits
{
    template<
        typename T_Species,
        typename T_UnspecifiedSpecies
    >
    struct SpeciesEligibleForSolver<
        T_Species,
        plugins::xrayDiffraction::XrayDiffraction< T_UnspecifiedSpecies >
    >
    {
        using Frame = typename T_Species::FrameType;

        using RequiredIdentifiers = MakeSeq_t<
            localCellIdx,
            position<>,
            weighting
        >;

        using type = typename pmacc::traits::HasIdentifiers<
            Frame,
            RequiredIdentifiers
        >::type;
    };

} // namespace traits
} // namespace particles
} // namespace picongpu
