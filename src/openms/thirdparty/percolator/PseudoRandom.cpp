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

#include "PseudoRandom.h"

namespace OpenMS { namespace Internal { namespace Percolator {

uint64_t PseudoRandom::seed_ = 1u;

// Park–Miller random number generator
// from wikipedia
unsigned long PseudoRandom::lcg_rand() {
  seed_ = (seed_ * 279470273u) % 4294967291u;
  return seed_;
}

// Generates a random double between 0 and 1
double PseudoRandom::lcg_uniform_rand() {
  return (double)PseudoRandom::lcg_rand() / ((double)PseudoRandom::kRandMax + (double)1);
}

// OpenMS: the state after `draws` calls of lcg_rand(state). The first draw is made as lcg_rand
// makes it (a state set by setSeed may exceed the modulus); from then on the state is below the
// modulus, and each draw multiplies it by 279470273 modulo 4294967291 (no product of two numbers
// below 2^32 overflows 64 bits), so the remaining draws multiply it by 279470273^(draws - 1).
uint64_t PseudoRandom::advance(uint64_t state, uint64_t draws) {
  if (draws == 0) return state;
  PseudoRandom::lcg_rand(state);
  uint64_t factor = 1u, power = 279470273u;
  for (uint64_t k = draws - 1; k != 0; k >>= 1) {
    if (k & 1u) factor = (factor * power) % 4294967291u;
    power = (power * power) % 4294967291u;
  }
  return (state * factor) % 4294967291u;
}

}}}  // namespace OpenMS::Internal::Percolator
