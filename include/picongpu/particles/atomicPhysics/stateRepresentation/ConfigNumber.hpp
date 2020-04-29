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
 * This file defines the atomic state representation by numbering Hydrogen-like
 * superconfigurations.
 *
 * The Hydrogen-like superconfiguration here is specified by the occupation
 * vector N-vector, a listing of every level n's occupation number N_n, with n:
 * 1 < n < n_max.
 *
 * These different superconfigurations are organized in a combinatorial table
 * and numbered starting with the completly ionized state.
 *
 * # |N_1, |N_2, |N_3, ...
 * 0 |0    |0    |0
 * 1 |1    |0    |0
 * 2 |2    |0    |0
 * 3 |0    |1    |0
 * 4 |1    |1    |0
 * 5 |2    |1    |0
 * 6 |0    |2    |0
 * 7 |1    |2    |0
 * 8 |2    |2    |0
 * 9 |0    |3    |0
 * ...
 *
 * analytical formula:
 * # ... configNumber assigned
 * # = N_1 *(g(0) + 1)+ N_2 * (g(1)+1) + N_3 * (g(2)+1) * (g(1) + 1)) + ...
 * # = Sum_{n=1}^{n_max}[ N_n * Produkt_{i=1}^{n} (g(i-1) + 1) ]
 *
 * g(n) ... maximum number of electrons in a given level n
 * g(n) = 2 * n^2
 * quick reference:
 * https://en.wikipedia.org/wiki/Electron_shell#Number_of_electrons_in_each_shell
 *
 * Note: a superconfiguration only stores occupation numbers, NOT spin or
 *  angular momentumm, due to memory constrains.
 *
 * further information see:
 *  https://en.wikipedia.org/wiki/Hydrogen-like_atom#Schr%C3%B6dinger_solution
 *
 * @todo 2020-07-01 BrianMarre: implement usage of current charge to account for
 * one more level than actually saved, => n_max effective = n_max + 1
 */

#pragma once


#include <pmacc/math/Vector.hpp>
#include <pmacc/algorithms/math/defines/comparison.hpp>

namespace picongpu
{
namespace particles
{
namespace atomicPhysics
{
namespace stateRepresentation
{

template< typename T_DataType, uint8_t T_NumberLevels, uint8_t T_ChargeNumber>
class ConfigNumber
{
 /* this class implements the actual storage of the configuration
 *
 * T_Numberlevels ... n_max
 * for convenience of usage and modularity, methods to convert the configNumber
 * to a occupation vector and convert a occuptation vector to the corresponding
 * configuration number are implemented.
 *
 */
    T_DataType configNumber;    // storage of actual configNumber

    uint16_t g(uint8_t n)
    {
        // cast necessary to prevent overflow in n^2 calculation
        return (static_cast<uint16_t>(n) * static_cast<uint16_t>(n) * 2);
    }

    T_DataType stepLength(uint8_t n)
    {
        T_DataType result = 1;

        for (uint8_t i = 1u; i < n; i++)
        {
            result *= pmacc::algorithms::math::min(this->g(i), T_ChargeNumber) + 1;
        }

        return result;
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
        PMACC_ASSERT_MSG(
            N < stepLength(T_NumberLevels + 1),
            "configurationNumber N larger than largest possible for"
            " T_NumberLevels"
        );

        this->configNumber = N;
    }

    ConfigNumber(
        pmacc::math::Vector< uint16_t, T_NumberLevels > levelVector
        )
    {
    /** constructor using a given occupation number vector to initialise.
    *
    * Uses the formula in file descripton. Assumes index of vector corresponds
    * to n-1, n ... principal quantum number.
    */

        /* stepLength ... number of table entries per occupation number VALUE of
        * the current principal quantum number n.
        */
        T_DataType stepLength = 1;

        this->configNumber = 0;

        for(uint8_t n=0u; n < T_NumberLevels; n++)
        {
            PMACC_ASSERT_MSG(
                this->g(n+1) >= levelVector[n],
                "occuation number too large"
            );

            /* BEWARE: n here equals n-1 in formula in file documentation,
            *
            * since for loop starts with n=0 instead of n=1,
            */
            stepLength *= std::min(this->g(n) + 1);
            this->configNumber += *levelVector[n] * stepLength;
        }
    }

    uint8_t numberLevels()
    {
    /** returns number of levels, n_max, used for configNumber
    */
        return T_NumberLevels;
    }

    operator pmacc::math::Vector< uint16_t, T_NumberLevels >()
    {
    /** B() operator converts configNumber B into an occupation number vector
    *
    * Index of result vector corresponds to principal quantum number n -1.
    *
    * exploits that for largest whole number k for a given configNumber N,
    * such that k * stepLength <= N, k is equal to the occupation number of
    * n_max.
    *
    * stepLength ... number of table entries per occupation number VALUE of
    * the current principal quantum number n.
    *
    * This is used recursively to determine all occupation numbers.
    * further information: see master thesis of Brian Marre
    */
        pmacc::math::Vector< uint16_t, T_NumberLevels > result =
            pmacc::math::Vector<uint16_t, T_NumberLevels>::create( 0 );

        T_DataType stepLength;
        T_DataType N;

        N = this->configNumber;

        // BEWARE: for loop counts down, strating with n_max
        for (uint8_t n = T_NumberLevels; n >= 1; n--)
        {
            // calculate current stepLength
            stepLength = this->stepLength(n);

            // get occupation number N_n by getting largest whole number factor
            *result[n-1] = static_cast<uint16_t>( N / stepLength );

            // remove contribution of current N_n
            N -= stepLength * (*result[n-1]);
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

// defines what datatype is to be used to save the data in this object
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

// defines how many independent components are saved in the object
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