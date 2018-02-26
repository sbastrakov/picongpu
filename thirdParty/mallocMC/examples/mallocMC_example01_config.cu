/*
  mallocMC: Memory Allocator for Many Core Architectures.
  https://www.hzdr.de/crp

  Copyright 2014 Institute of Radiation Physics,
                 Helmholtz-Zentrum Dresden - Rossendorf

  Author(s):  Carlchristian Eckert - c.eckert ( at ) hzdr.de

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/

#pragma once

#include <boost/mp11/integral.hpp>

// basic files for mallocMC
#include "src/include/mallocMC/mallocMC_hostclass.hpp"

// Load all available policies for mallocMC
#include "src/include/mallocMC/CreationPolicies.hpp"
#include "src/include/mallocMC/DistributionPolicies.hpp"
#include "src/include/mallocMC/OOMPolicies.hpp"
#include "src/include/mallocMC/ReservePoolPolicies.hpp"
#include "src/include/mallocMC/AlignmentPolicies.hpp"
    


// configurate the CreationPolicy "Scatter" to modify the default behaviour
struct ScatterHeapConfig : mallocMC::CreationPolicies::Scatter<>::HeapProperties{
    typedef boost::mp11::mp_int<4096> pagesize;
    typedef boost::mp11::mp_int<8>    accessblocks;
    typedef boost::mp11::mp_int<16>   regionsize;
    typedef boost::mp11::mp_int<2>    wastefactor;
    typedef boost::mp11::mp_false     resetfreedpages;
};

struct ScatterHashConfig : mallocMC::CreationPolicies::Scatter<>::HashingProperties{
    typedef boost::mp11::mp_int<38183> hashingK;
    typedef boost::mp11::mp_int<17497> hashingDistMP;
    typedef boost::mp11::mp_int<1>     hashingDistWP;
    typedef boost::mp11::mp_int<1>     hashingDistWPRel;
};

// configure the DistributionPolicy "XMallocSIMD"
struct XMallocConfig : mallocMC::DistributionPolicies::XMallocSIMD<>::Properties {
  typedef ScatterHeapConfig::pagesize pagesize;
};

// configure the AlignmentPolicy "Shrink"
struct ShrinkConfig : mallocMC::AlignmentPolicies::Shrink<>::Properties {
  typedef boost::mp11::mp_int<16> dataAlignment;
};

// Define a new allocator and call it ScatterAllocator
// which resembles the behaviour of ScatterAlloc
typedef mallocMC::Allocator< 
  mallocMC::CreationPolicies::Scatter<ScatterHeapConfig, ScatterHashConfig>,
  mallocMC::DistributionPolicies::XMallocSIMD<XMallocConfig>,
  mallocMC::OOMPolicies::ReturnNull,
  mallocMC::ReservePoolPolicies::SimpleCudaMalloc,
  mallocMC::AlignmentPolicies::Shrink<ShrinkConfig>
  > ScatterAllocator;
