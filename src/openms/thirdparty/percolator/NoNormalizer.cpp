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
#include <iostream>
#ifdef WIN32
#include <float.h>
#define isfinite _finite
#endif
#include <cmath>
#include <vector>
#include <string>

using namespace std;

#include "Normalizer.h"
#include "NoNormalizer.h"


NoNormalizer::NoNormalizer() {
}

NoNormalizer::~NoNormalizer() {
}

void NoNormalizer::unnormalizeweight(const std::vector<double>& in, 
    std::vector<double>& out) const {
  unsigned int i = 0;
  for (; i < numFeatures+1; i++) {
    out[i] = in[i];
  }
}

void NoNormalizer::normalizeweight(const std::vector<double>& in, 
    std::vector<double>& out) const {
  size_t i = 0;
  for (; i < numFeatures+1; i++) {
    out[i] = in[i];
  }
}

void NoNormalizer::setSet(std::vector<double*>& /* featuresV */,
                            std::vector<double*>& /* rtFeaturesV */, size_t nf,
                            size_t nrf) {
  numFeatures = nf;
  numRetentionFeatures = nrf;
  sub.resize(nf + nrf, 0.0);
  div.resize(nf + nrf, 0.0);
  size_t ix;
  for (ix = 0; ix < numFeatures; ++ix) {
    sub[ix] = 0.;
    div[ix] = 1.0;
  }
}


void NoNormalizer::updateSet(vector<double*> & /* featuresV */, size_t /* offset */,
                               size_t numFeatures) {
  size_t ix;
  for (ix = 0; ix < numFeatures; ++ix) {
    sub[ix] = 0.;
    div[ix] = 1.0;
  }
}
