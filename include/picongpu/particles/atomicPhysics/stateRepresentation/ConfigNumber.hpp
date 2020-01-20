/* Copyright 2019-2020 Brian Marre
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

/** @file
 *
 * This file defines the atomic state representation by numbering Hydrogen
 * superconfigurations.
 *
 * The Hydrogen superconfiguration is specified by the occupation vector
 * N-vector, a listing of every level n's occupation number N_n,
 * Hydrogen 1 < n < n_max.
 * These different configs. are organized in a combinatorial table and numbered
 * starting with the completly ionized state.
 *
 * # |N_1, |N_2, |N_3; ...
 * 0 |0    |0    |0
 * 1 |1    |0    |0
 * 2 |2    |0    |0
 * 3 |0    |1    |0
 * 4 |1    |1    |0
 * 5 |2    |1    |0
 * 6 |0    |2    |0
 *
 * #...config Number assigned,
 * g(n)= 2 * n^2 ...maximum number of electrons in a given level n
 * # = N_1 *(0 + 1)+ N_2 * (g(1)+1) + N_3 * (g(2)+1) * (g(1) + 1)) + ...
 *
 * # = Sum_n=1^n_max( N_n * Produkt_i=1^n (g(i-1) +1) )
 *
 *@todo 2020-07-01 BrianMarre: implement usage of current charge to account for
 * one more than actually saved levels, => n_max effective = n_max + 1
 */

#pragma once


#include <pmacc/math/Vector.hpp>

namespace picongpu
{
namespace particles
{
namespace atomicPhysics
{
namespace stateRepresentation
{

template< typename T_DataType, uint8_t T_NumberLevels >
class ConfigNumber
{
 /* this class implements the actual storage of the Configuration
 *
 *T_Numberlevels ... n_max
 *for convenience of usage and modularity methods to convert the Config Number
 *to a occupation vector and convert a occuptation vector to the corresponding
 *configuration number are implemented.
 *
 */
    T_DataType configNumber;    //storage of actual ConfigNumber

    uint16_t g(uint8_t n)
    {
        return static_cast<unsigned short int>(n) * n * 2;
    }

public:
    ConfigNumber(
        T_DataType N = 0u
        )
    {
        PMACC_ASSERT_MSG(
            N >= 0,
            "negative configurationNumbers are not defined"
        );

        this->configNumber = N;
    }
    ConfigNumber(
        pmacc::math::Vector< uint8_t, T_NumberLevels > levelVector
        )
    {
        T_DataType stepLength = 1;
        this->configNumber = 0;

        for(uint8_t n=0u; n < T_NumberLevels; n++) {
            PMACC_ASSERT_MSG(
                this->g(n) >= *levelVector[n],
                "occuation number too large"
            );

            stepLength *= this->g(n) + 1;
            this->configNumber += *levelVector[n] * stepLength;
        }
    }

    uint8_t numberLevels()
    {
        return T_NumberLevels;
    }

    operator pmacc::math::Vector< uint8_t, T_NumberLevels >()
    {
        pmacc::math::Vector< uint8_t, T_NumberLevels > result =
            pmacc::math::Vector<uint8_t, T_NumberLevels>::create( 0 );

        T_DataType product;
        T_DataType N;

        N = this->configNumber;

        for (uint8_t n = T_NumberLevels; n >= 1; n--)
        {
            product = 1;
            for (uint8_t i = 1u; i < n; i++)
            {
                product *= ( this->g(i) + 1 );
            }

            *result[n-1] = static_cast<uint8_t>( N / product );

            N -= product * (*result[n-1]);
        }

        return result;
    }
};

} // namespace stateRepresentation
} // namespace atomicPhysics
} // namespace particles
} // namespace picongpu

/** this specfies how an object of the ConfigNumber class can be written to an
 * external file for storage.
*/

namespace pmacc
{
namespace traits
{

//defines what datatype is to be used to save the data in this object
template< typename T_DataType, uint8_t T_NumberLevels >
struct GetComponentsType<
    picongpu::particles::atomicPhysics::stateRepresentation::ConfigNumber<
        T_DataType,
        T_NumberLevels
    >,
    false
>
{
    using type = T_DataType;
};

//defines how many independent components are saved in the object
template< typename T_DataType, uint8_t T_NumberLevels >
struct GetNComponents<
    picongpu::particles::atomicPhysics::stateRepresentation::ConfigNumber<
        T_DataType,
        T_NumberLevels
    >,
    false
>
{
    static constexpr uint32_t value = 1u;
};

} // namespace traits
} // namespace pmacc
