/* Copyright 2014-2020 Axel Huebl, Felix Schmitt, Heiko Burau, Rene Widera
 *                     Benjamin Worpitz, Sergei Bastrakov
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

#include <pmacc/types.hpp>
#include "picongpu/simulation_defines.hpp"
#include "picongpu/plugins/adios/ADIOSWriter.def"
#include "picongpu/plugins/misc/ComponentNames.hpp"
#include "picongpu/traits/IsFieldDomainBound.hpp"

#include <pmacc/particles/frame_types.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/dimensions/DataSpace.hpp>
#include <pmacc/dimensions/GridLayout.hpp>
#include "picongpu/simulation/control/MovingWindow.hpp"

#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>

#include <string>
#include <sstream>
#include <stdexcept>


namespace picongpu
{

namespace adios
{

/**
 * Helper class for ADIOS plugin to load fields from parallel ADIOS BP files.
 */
class RestartFieldLoader
{
public:
    template<class Data>
    static void loadField(
        Data& field,
        const uint32_t numComponents,
        std::string objectName,
        ThreadParams *params,
        const bool isDomainBound
    )
    {
        log<picLog::INPUT_OUTPUT > ("Begin loading field '%1%'") % objectName;

        const auto componentNames = plugins::misc::getComponentNames( numComponents );
        const DataSpace<simDim> field_guard = field.getGridLayout().getGuard();

        const pmacc::Selection<simDim>& localDomain = Environment<simDim>::get().SubGrid().getLocalDomain();

        using ValueType = typename Data::ValueType;
        field.getHostBuffer().setValue(ValueType::create(0.0));

        DataSpace<simDim> domain_offset = localDomain.offset;
        DataSpace<simDim> local_domain_size = params->window.localDimensions.size;
        int elementCount = params->window.localDimensions.size.productOfComponents();
        bool useLinearIdxAsDestination = false;

        /* Patch for non-domain-bound fields
        * This is an ugly fix to allow output of reduced 1d PML buffers,
        * that are the same size on each domain.
        * This code is to be replaced with the openPMD output plugin soon.
        */
        if( !isDomainBound )
        {
            auto const field_layout = params->gridLayout;
            auto const field_no_guard = field_layout.getDataSpaceWithoutGuarding();
            elementCount = field_no_guard.productOfComponents();
            auto const & gridController = Environment<simDim>::get().GridController();
            auto const rank = gridController.getGlobalRank();
            domain_offset = DataSpace<simDim>::create( 0 );
            domain_offset[ 0 ] = rank * elementCount;
            local_domain_size = DataSpace<simDim>::create( 1 );
            local_domain_size[ 0 ] = elementCount;
            useLinearIdxAsDestination = true;
            
            std::cout << "\n\nPML start\nelementCount = " << elementCount << "\n";
        }

        auto destBox = field.getHostBuffer().getDataBox();
        for (uint32_t n = 0; n < numComponents; ++n)
        {
            // Read the subdomain which belongs to our mpi position.
            // The total grid size must match the grid size of the stored data.
            log<picLog::INPUT_OUTPUT > ("ADIOS: Read from domain: offset=%1% size=%2%") %
                domain_offset % local_domain_size;

            std::stringstream datasetName;
            datasetName << params->adiosBasePath << ADIOS_PATH_FIELDS << objectName;
            if (numComponents > 1)
                datasetName << "/" << componentNames[n];

            log<picLog::INPUT_OUTPUT > ("ADIOS: Read from field '%1%'") %
                datasetName.str();

            ADIOS_VARINFO* varInfo = adios_inq_var( params->fp, datasetName.str().c_str() );
            if( varInfo == nullptr )
            {
                std::string errMsg( adios_errmsg() );
                if( errMsg.empty() ) errMsg = '\n';
                std::stringstream s;
                s << "ADIOS: error at adios_inq_var '"
                  << "' (" << adios_errno << ") in "
                  << __FILE__ << ":" << __LINE__ << " " << errMsg;
                throw std::runtime_error(s.str());
            }
            uint64_t start[varInfo->ndim];
            uint64_t count[varInfo->ndim];
            for(int d = 0; d < varInfo->ndim; ++d)
            {
                /* \see adios_define_var: z,y,x in C-order */
                start[d] = domain_offset.revert()[d];
                count[d] = local_domain_size.revert()[d];
            }

            std::cout << "start = " << start[0] << " " << start[1] << " " << start[2] << "\n";
            std::cout << "count = " << count[0] << " " << count[1] << " " << count[2] << "\n";
            std::cout << "local_domain_size = " << local_domain_size[0] << " " << local_domain_size[1] << " " << local_domain_size[2] << "\n";
            
            ADIOS_SELECTION* fSel = adios_selection_boundingbox( varInfo->ndim, start, count );

            /* specify what we want to read, but start reading at below at
             * `adios_perform_reads` */
            log<picLog::INPUT_OUTPUT > ("ADIOS: Allocate %1% elements") %
                local_domain_size.productOfComponents();

            /// \todo float_X should be some kind of gridBuffer's GetComponentsType<ValueType>::type
            float_X* field_container = new float_X[local_domain_size.productOfComponents()];
            /* magic parameters (0, 1): `from_step` (not used in streams), `nsteps` to read (must be 1 for stream) */
            log<picLog::INPUT_OUTPUT > ("ADIOS: Schedule read from field (%1%, %2%, %3%, %4%)") %
                                        params->fp % fSel % datasetName.str() % (void*)field_container;

            // avoid deadlock between not finished pmacc tasks and mpi calls in adios
            __getTransactionEvent().waitForFinished();
            ADIOS_CMD(adios_schedule_read( params->fp, fSel, datasetName.str().c_str(), 0, 1, (void*)field_container ));

            /* start a blocking read of all scheduled variables */
            ADIOS_CMD(adios_perform_reads( params->fp, 1 ));

            ///int elementCount = params->window.localDimensions.size.productOfComponents();

            #pragma omp parallel for
            for (int linearId = 0; linearId < elementCount; ++linearId)
            {
                DataSpace<simDim> destIdx;
                if( useLinearIdxAsDestination )
                {
                    destIdx[ 0 ] = linearId;
                }
                else
                {
                    /* calculate index inside the moving window domain which is located on the local grid*/
                    destIdx = DataSpaceOperations<simDim>::map(params->window.localDimensions.size, linearId);
                    /* jump over guard and local sliding window offset*/
                    destIdx += field_guard + params->localWindowToDomainOffset;
                }
                destBox(destIdx)[n] = field_container[linearId];
            }

            __deleteArray(field_container);
            adios_selection_delete(fSel);
            adios_free_varinfo(varInfo);
        }

        field.hostToDevice();

        __getTransactionEvent().waitForFinished();

        log<picLog::INPUT_OUTPUT > ("ADIOS: Read from domain: offset=%1% size=%2%") %
            domain_offset % local_domain_size;
        log<picLog::INPUT_OUTPUT > ("ADIOS: Finished loading field '%1%'") % objectName;
    }

};

/**
 * Helper class for ADIOSWriter (forEach operator) to load a field from ADIOS
 *
 * @tparam T_Field field class to load
 */
template< typename T_Field >
struct LoadFields
{
public:

    HDINLINE void operator()(ThreadParams* params)
    {
#ifndef __CUDA_ARCH__
        DataConnector &dc = Environment<>::get().DataConnector();
        ThreadParams *tp = params;

        /* load field without copying data to host */
        auto field = dc.get< T_Field >( T_Field::getName(), true );
        tp->gridLayout = field->getGridLayout();

        /* load from ADIOS */
        bool const isDomainBound = traits::IsFieldDomainBound< T_Field >::value;
        RestartFieldLoader::loadField(
            field->getGridBuffer(),
            (uint32_t)T_Field::numComponents,
            T_Field::getName(),
            tp,
            isDomainBound
        );

        dc.releaseData(T_Field::getName());
#endif
    }

};

using namespace pmacc;

} /* namespace adios */
} /* namespace picongpu */
