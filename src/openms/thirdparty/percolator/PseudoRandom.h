/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef PSEUDO_RANDOM_H
#define PSEUDO_RANDOM_H


#include <stdint.h>

namespace OpenMS { namespace Internal { namespace Percolator {

/*
* Random is a helper class generating pseudo random numbers starting from a seed
*
* Here are some usefull abbreviations:
* LCG - Linear Congruential Generator
*
*/
class PseudoRandom {
 public:
  inline static void setSeed(unsigned long s) { seed_ = s; }
  static unsigned long lcg_rand();
  static double lcg_uniform_rand();
  const static uint64_t kRandMax = 4294967291u;
  // OpenMS: the numbers of lcg_rand and lcg_uniform_rand from a state of the caller's (to draw
  // parts of the sequence in parallel), the generator's state, and the state after a number of
  // draws (jumps ahead without drawing)
  inline static unsigned long lcg_rand(uint64_t& state) {
    state = (state * 279470273u) % 4294967291u;
    return state;
  }
  inline static double lcg_uniform_rand(uint64_t& state) {
    return (double)PseudoRandom::lcg_rand(state) / ((double)PseudoRandom::kRandMax + (double)1);
  }
  inline static uint64_t getSeed() { return seed_; }
  static uint64_t advance(uint64_t state, uint64_t draws);
 protected:
  static uint64_t seed_;
};


}}}  // namespace OpenMS::Internal::Percolator

#endif // RANDOM_H
