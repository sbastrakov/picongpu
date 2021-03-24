/* Copyright 2020-2021 Brian Marre
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
 */

#include "picongpu/simulation_defines.hpp"
#include "picongpu/particles/atomicPhysics/AtomicRate.hpp"

#include <utility>
#include <pmacc/algorithms/math.hpp>

// debug only
#include <iostream>
#include <cmath>

#pragma once

namespace picongpu
{
    namespace particles
    {
        namespace atomicPhysics
        {
            namespace electronDistribution
            {
                namespace histogram2
                {
                    template<
                        uint32_t T_minOrderApprox, // half of starting order of error approximation,
                                                   // in crosssection + velocity and 0th order in
                                                   // electron density
                        uint32_t T_maxOrderApprox, // half of maximum order of error approximation,
                                                   // T_minOrderapprox <= a <= T_maxOrderApprox,
                                                   //  2*a ... order of a term of the rate approximation
                        uint32_t T_numSamplePoints, // number of sample points used for numerical differentiation
                                                    // >= 2*a + 1
                        typename T_WeightingGen, // type of numerical differentiation method used
                        typename T_AtomicRate, // atomic rate functor
                        typename T_AtomicDataBox,
                        typename T_ConfigNumberDataType>
                    class RateRelativeError
                    {
                        // TODO: add pmacc assert to check T_numSamplePoints >= 2*T_maxOrderApprox + 1
                        using AtomicRate = T_AtomicRate;
                        // necessary for velocity derivation

                        // storage of weights for later use
                        //float_X weights[T_numSamplePoints * (2u * T_maxOrderApprox + 1u)];
                        float_X weights[ 5u*5u ];

                        // return unitless
                        /** returns the k-th of T_numNodes chebyshev nodes x_k, interval [-1, 1]
                         *
                         * x_k = cos( (2k-1)/(2n) * pi )
                         *
                         * @tparam T_numNodes ... how many chebyshev nodes are requested
                         * @param k ... index of chebyshev Nodes
                         *
                         * BEWARE: k=0 is NOT allowed, k \in {1, ... , T_numSamplePoints}
                         * BEWARE: max({ x_k }) = x_1 and min({ x_k }) = x_(T_numSamplePoints)
                         *
                         * see https://en.wikipedia.org/wiki/Chebyshev_nodes for more information
                         */
                        template<uint32_t T_numNodes>
                        DINLINE static float_X chebyshevNodes(uint32_t const k)
                        {
                            // check for bounds on k
                            // return boundaries if outside of boundary
                            if(k < 1u)
                                return -1._X;
                            if(k > T_numNodes)
                                return 1._X;

                            return math::cos<float_X>(
                                static_cast<float_X>(2u * k - 1u) / static_cast<float_X>(2 * T_numNodes)
                                * picongpu::PI);
                        }


                        // return unitless
                        DINLINE static float_64 fak(const uint32_t n)
                        {
                            float_64 result = 1u;

                            for(uint32_t i = 1u; i <= n; i++)
                            {
                                result *= static_cast<float_64>(i);
                            }

                            return result;
                        }

                       DINLINE static float_X samplePoints( uint8_t j)
                       {
                            if(j < 4)
                            {
                                return -1._X/2._X + float_X(j);
                            }
                            return 0._X;
                       }

                    public:
                        // calculates weightings once
                        DINLINE void init()
                        {
                            //float_X samplePoints[T_numSamplePoints];

                            /*
                            // include reference point
                            samplePoints[0] = 0._X;

                            for(uint32_t i = 1u; i < T_numSamplePoints; i++)
                            {
                                samplePoints[i] = chebyshevNodes<(T_numSamplePoints - 1u)>(i);
                            }

                            // calculate weightings for each order derivative until maximum
                            for(uint32_t i = 0u; i <= (2u * T_maxOrderApprox); i++) // i ... order of derivation
                            {
                                for(uint32_t j = 0u; j < T_numSamplePoints; j++) // j ... sample point
                                {
                                    this->weights[j + i * T_numSamplePoints] = T_WeightingGen::template weighting<
                                        T_numSamplePoints,
                                        T_numSamplePoints / 2u + 2u>(i, j, samplePoints);
                                }
                            }*/
                            //T_WeightingGen::template test();

                            // debug workaround
                            // 0th order derivative
                            this->weights[0u] = float_X(35./128.);
                            this->weights[1u] = float_X(35./32.);
                            this->weights[2u] = float_X(-35./64.);
                            this->weights[3u] = float_X(7./32.);
                            this->weights[4u] = float_X(-5./128.);

                            // 1th order derivative
                            this->weights[0u+1u*5u] = float_X(-11./12.);
                            this->weights[1u+1u*5u] = float_X(17./24.);
                            this->weights[2u+1u*5u] = float_X(3./8.);
                            this->weights[3u+1u*5u] = float_X(-5./24.);
                            this->weights[4u+1u*5u] = float_X(1./24.);

                            // 2th order derivative
                            this->weights[0u+2u*5u] = float_X(43./24.);
                            this->weights[1u+2u*5u] = float_X(-14./3.);
                            this->weights[2u+2u*5u] = float_X(17./4.);
                            this->weights[3u+2u*5u] = float_X(-5./3.);
                            this->weights[4u+2u*5u] = float_X(7./24.);

                            // 3th order derivative
                            this->weights[0u+3u*5u] = float_X(-2.);
                            this->weights[1u+3u*5u] = float_X(7.);
                            this->weights[2u+3u*5u] = float_X(-9.);
                            this->weights[3u+3u*5u] = float_X(5.);
                            this->weights[4u+3u*5u] = float_X(-1.);

                            // 4th order derivative
                            this->weights[0u+4u*5u] = float_X(1);
                            this->weights[1u+4u*5u] = float_X(-4);
                            this->weights[2u+4u*5u] = float_X(6);
                            this->weights[3u+4u*5u] = float_X(-4);
                            this->weights[4u+4u*5u] = float_X(1);

                            /*
                            // debug only
                            printf("samplePoints {%f, %f, %f}\n",
                                samplePoints[0],
                                samplePoints[1],
                                samplePoints[2]
                                );*/
                        }


                        // return unit: [1/s]/[ 1/(m^3 * J) ], error per electron density
                        template<typename T_Acc>
                        DINLINE float_X operator()(
                            T_Acc& acc,
                            float_X dE, // unit: ATOMIC_UNIT_ENERGY
                            float_X E, // unit: ATOMIC_UNIT_ENERGY
                            T_AtomicDataBox atomicDataBox,
                            bool debug=true) const
                        {
                            float_X result = 0._X; // unit: (1/s) / (m^3/J)

                            // sanity checks
                            // dE > 0, otherwise rate relative error not defined (dE==0),
                            // or incorrectly calculated (dE < 0)
                            if(dE <= 0._X)
                            {
                                printf("WARNING: rate relative error not defined, for binWidth: %f\n", dE);
                            }

                            // E > 0, required by pyhsics due to definition of energy
                            if(E <= 0._X)
                            {
                                printf("WARNING: rate relative error not defined, for central bin Energy: %f\n", E);
                            }

                            // debug only
                            uint16_t loopCounter = 0u;

                            // a ... order of rate approximation
                            for(uint32_t a = T_minOrderApprox; a <= T_maxOrderApprox; a++)
                            {
                                // o ... derivative order of crossection sigma
                                for(uint32_t o = 0u; o <= 2 * a; o++)
                                {
                                    // debug only
                                    loopCounter++;

                                    /*std::cout << "loopCounter " << loopCounter;
                                    std::cout << ", fak(" << o <<") " << fak(o);

                                    std::cout << ", fak(" << (2u * a - o) << ") " << fak(2u * a - o);
                                    std::cout << ", sigmaDerivative(E:" << E <<", dE:" << dE << ", o:" << o <<") " <<
                                        sigmaDerivative(acc, E, dE, o, atomicDataBox);
                                    std::cout << ", velocityDerivative(E:" << E << ", dE:" << dE << ", k:"<< (2u * a - o) << ") "<<
                                        velocityDerivative(E, dE, 2u * a - o) * 1._X / (2u * a + 1u);
                                    std::cout << " conversion factors " << math::pow(
                                                  dE * picongpu::SI::ATOMIC_UNIT_ENERGY / 2._X,
                                                  static_cast<float_X>(2u * a + 1u));
                                    std::cout << std::endl;
                                    // printf("loopCounter %i, sigmaDerivative %f\n", loopCounter, sigmaDerivative(acc,
                                    // E, dE, o, atomicDataBox));*/

                                    // taylor expansion of rate integral, see my master thesis
                                    // \sum{ 1/(o! (2a-o)!) * sigma'(o) * v'(2a-o) * 1/(2a+1) * (dE/2)^(2a+1) * 2}
                                    // 1/unitless * m^2/J^o * m/s * 1/J^(2*a-o) * unitless * J^(2*a+1) * unitless
                                    // = m^3/J^(2a) * 1/s * J^(2a+1) = J * m^3 * 1/s
                                    result += 1.0_X / (fak(o) * fak(2u * a - o))
                                        * sigmaDerivative(acc, E, dE, o, atomicDataBox, debug)
                                        * velocityDerivative(E, dE, 2u * a - o) * 1._X / (2u * a + 1u)
                                        * math::pow(
                                                  dE * picongpu::SI::ATOMIC_UNIT_ENERGY / 2._X,
                                                  static_cast<float_X>(2u * a + 1u))
                                        * 2._X; // unit: J * m^3 * 1/s
                                }
                            }
                            // debug only
                            // printf("    deltaE: %f\n", dE);
                            // printf("    relative error: %f\n", result);
                            return result;
                        }

                    private:
                        // @param energy ... unit: atomic energy units
                        // return unit: m/s, SI
                        // TODO: remove assumption of electron mass
                        DINLINE static float_X velocity(float_X energy // unit: atomic energy unit
                        )
                        {
                            constexpr float_X mass_e_SI = picongpu::SI::ELECTRON_MASS_SI;
                            constexpr float_X c_SI = picongpu::SI::SPEED_OF_LIGHT_SI;

                            float_X restEnergy_SI = mass_e_SI * math::pow(c_SI, 2.0_X);

                            // v = sqrt( 1 - (m^2*c^4)/(E^2) )
                            return picongpu::SI::SPEED_OF_LIGHT_SI
                                * math::sqrt(
                                       1._X
                                       - math::pow(
                                           1._X / (1._X + (energy * picongpu::SI::ATOMIC_UNIT_ENERGY) / restEnergy_SI),
                                           2.0_X));
                        }

                        // return unit: m/s * 1/J^m), SI
                        DINLINE float_X velocityDerivative(
                            float_X E, // unit: ATOMIC_UNIT_ENERGY
                            float_X dE, // unit: ATOMIC_UNIT_ENERGY
                            uint32_t m // order of derivative
                        ) const
                        {

                            float_X samplePoint; // uint: ATOMIC_UNIT_ENERGY
                            float_X weight; // unit: unitless
                            float_X velocityValue; // unit: m/s, SI

                            float_X result = 0._X; // unit: m/s, SI

                            // debug only
                            uint16_t loopCounter = 0u;

                            for(uint32_t j = 0u; j < /*T_numSamplePoints*/5u; j++) // j index of sample point
                            {
                                // debug only
                                loopCounter++;

                                samplePoint = E + samplePoints(j) * dE/2._X;
                                /*if( j == 0 )
                                {
                                    // first samplePoint, is by construction always = 0:
                                    samplePoint = E;
                                }
                                else
                                {
                                    // calculate chebyshev sample point
                                    samplePoint = E + chebyshevNodes<T_numSamplePoints -1u>(j) * dE / 2._X;
                                }*/

                                // get Weight from pregenerated table
                                weight = this->weights[m * /*T_numSamplePoints*/5u + j]; // unit: unitless

                                // velocity( [ ATOMIC_UNIT_ENERGY ] ) -> m/s, SI
                                velocityValue = velocity(samplePoint); // unit: m/s, SI

                                result += weight * velocityValue; // unit: m/s, SI

                                // debug only
                                std::cout << "loopCounter " << loopCounter << " samplePoint(m,j) " << samplePoint <<
                                    " weight " << weight << " velocityValue " << velocityValue << " result " <<
                                    result << std::endl;
                            }

                            // get correct return unit
                            // pow( 1/[ ATOMIC_UNIT_ENERGY * AU_To_J ] = 1/J, m )-> 1/(J^m)
                            result /= math::pow(
                                dE / 2._X * picongpu::SI::ATOMIC_UNIT_ENERGY,
                                static_cast<float_X>(m)); // unit: m/s * 1/J^m, SI

                            return result;
                        }


                        // return unit: m^2/J^m, SI
                        template<typename T_Acc>
                        DINLINE float_X sigmaDerivative(
                            T_Acc& acc,
                            float_X E, // central energy of bin, unit: ATOMIC_UNIT_ENERGY
                            float_X dE, // delta energy, unit: ATOMIC_UNIT_ENERGY
                            uint32_t o, // order of derivative, unitless
                            T_AtomicDataBox const atomicDataBox, // in file atomicData.hpp
                            bool debug=false
                        ) const
                        {
                            // debug only
                            // printf("E %f , dE %f, o %i\n", E, dE, o);

                            float_X samplePoint; // uint: ATOMIC_UNIT_ENERGY
                            float_X weight; // unit: unitless
                            float_X sigmaValue; // unit: m^2, SI

                            float_X result = 0._X; // unit: m^2, SI

                            // debug only
                            uint16_t loopCounter = 0u;

                            for(uint32_t j = 0u; j < /*T_numSamplePoints*/5u; j++) // j index of sample point
                            {
                                // debug only
                                loopCounter++;


                                samplePoint = E + samplePoints(j) * dE/2._X;
                                /*if( j == 0 )
                                {
                                    // first samplePoint, is by construction always = 0:
                                    samplePoint = E;
                                }
                                else
                                {
                                    // calculate chebyshev sample point
                                    samplePoint = E + chebyshevNodes<T_numSamplePoints -1u>(j) * dE / 2._X;
                                }*/

                                // get Weight from pregenerated table
                                // TODO: stuff this index arithmic into a method call/template
                                weight = this->weights[o * /*T_numSamplePoints*/5u + j]; // unit: unitless

                            // get crossection of the energy bin
                                float_X sigmaValue = AtomicRate::totalCrossSection(
                                    acc,
                                    samplePoint, // unit: ATOMIC_UNIT_ENERGY
                                    atomicDataBox,
                                    debug); // unit: m^2, SI

                                result += weight * sigmaValue; // unit: m^2, SI

                                // debug only
                                std::cout << "loopCounter " << loopCounter << " samplePoint(m,j) " << samplePoint <<
                                    " weight " << weight << " sigmaValue " << sigmaValue << " result " <<
                                    result << std::endl;
                            }

                            /*printf("step %i: weight %f, sigmaValue %f, result %f\n",
                                loopCounter, weight, sigmaValue, result);*/

                            // orignial version
                           /*sigmaValue = AtomicRate::totalCrossSection(
                              acc,
                              E
                                  + chebyshevNodes<numSamplePoints - 1u>(j + 1u) * dE
                                      / 2._X, // unit: ATOMIC_UNIT_ENERGY
                              atomicDataBox); // unit: m^2, SI*/


                            result /= math::pow(dE * picongpu::SI::ATOMIC_UNIT_ENERGY / 2._X, static_cast<float_X>(o));
                            //} m^2 / (J)^m

                            return result; // unit: m^2/J^m, SI
                        }
                    };

                } // namespace histogram2
            } // namespace electronDistribution
        } // namespace atomicPhysics
    } // namespace particles
} // namespace picongpu
