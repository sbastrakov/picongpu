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

#include "picongpu/plugins/xrayDiffraction/ComputeGlobalDomain.hpp"
#include "picongpu/plugins/xrayDiffraction/ReciprocalSpace.hpp"

#include <pmacc/Environment.hpp>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


namespace picongpu
{
namespace plugins
{
namespace xrayDiffraction
{
namespace detail
{

    /** Writer of X-ray diffraction results to files
     *
     * Intended to be instantiated once and used from a single rank in a
     * simulation
     */
    class Writer
    {
    public:

        /** Create a writer
         *
         * @param filePrefix output file prefix
         */
        Writer( std::string const & filePrefix );

        /** Write diffraction data to a file
         *
         * @param globalDomainResult result aggregated for the global domain
         * @param currentStep current time iteration
         */
        void write(
            ComputeGlobalDomain const & globalDomainResult,
            uint32_t currentStep
        ) const;

    private:

        //! Directory name inside the general output directory
        std::string const directoryName;

        //! Output file prefix
        std::string const filePrefix;

        /** Write diffraction intensity to a file
         *
         * @param diffractionIntensity diffraction intensity for all scattering vectors
         * @param reciprocalSpace reciprocal space
         * @param fileName file name
         */
        static void writeIntensity(
            std::vector< float_X > const & diffractionIntensity,
            ReciprocalSpace const & reciprocalSpace,
            std::string const & fileName
        );

        /** Write diffraction number of particles and macroparticles to a file
         *
         * @param numMacroparticles number of macroparticles
         * @param totalWeighting combined macroparticle weighting
         * @param fileName file name
         */
        static void writeLog( 
            int64_t numMacroparticles,
            float_64 totalWeighting,
            std::string const & fileName
        );

    };

    Writer::Writer( std::string const & filePrefix ):
        filePrefix( filePrefix ),
        directoryName( "xrayDiffraction" )
    {
        auto & fs = Environment< simDim >::get().Filesystem();
        fs.createDirectoryWithPermissions( directoryName );
    }

    void Writer::write(
        ComputeGlobalDomain const & globalDomainResult,
        uint32_t const currentStep
    ) const
    {
        auto const baseFileName = directoryName + filePrefix +
            "_" + std::to_string( currentStep );
        writeIntensity(
            globalDomainResult.diffractionIntensity,
            globalDomainResult.reciprocalSpace,
            baseFileName + ".dat"
        );
        writeLog(
            globalDomainResult.numMacroparticles,
            globalDomainResult.totalWeighting,
            baseFileName + ".log"
        );
    }

    void Writer::writeIntensity(
        std::vector< float_X > const & diffractionIntensity,
        ReciprocalSpace const & reciprocalSpace,
        std::string const & fileName
    )
    {
        auto ofile = std::ofstream{ fileName.c_str() };
        if( !ofile )
            std::cerr << "Could not open file [" << fileName
            << "] for output, disable plugin output.\n";
        else
        {
            ofile << diffractionIntensity.size() << "\n";
            ofile << "# qx qy qz intensity \n";
            for( size_t i = 0; i < diffractionIntensity.size(); i++ )
                ofile << reciprocalSpace.getValue( i ) << " "
                << diffractionIntensity[ i ] << "\n";
        }
    }

    void Writer::writeLog( 
        int64_t numMacroparticles,
        float_64 totalWeighting,
        std::string const & fileName
    )
    {
        auto ofile = std::ofstream{ fileName.c_str() };
        if( !ofile )
            std::cerr << "Can't open file [" << fileName
            << "] for X-ray diffraction output, disable plugin output.\n";
        else
        {
            ofile << "Number of macroparticles:"
                << " " << numMacroparticles << "\n";
            ofile << "Number of particles (total weighting):"
                << " " << totalWeighting << "\n";
        }
    }

} // namespace detail
} // namespace xrayDiffraction
} // namespace plugins
} // namespace picongpu
