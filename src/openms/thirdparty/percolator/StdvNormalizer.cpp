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
#include <cmath>
#include <limits>
#include <vector>
#include <string>


#include "Normalizer.h"
#include "StdvNormalizer.h"
#include "DataSet.h"
#include "FeatureNames.h"
#include "Globals.h"

namespace OpenMS { namespace Internal { namespace Percolator {

using namespace std;


namespace {

// The smallest and the largest value of each column. A column is constant
// exactly if the two are equal.
struct ColumnRanges {
  explicit ColumnRanges(size_t numColumns)
      : minV(numColumns, std::numeric_limits<double>::infinity()),
        maxV(numColumns, -std::numeric_limits<double>::infinity()) {}
  void add(size_t ix, double value) {
    if (value < minV[ix]) minV[ix] = value;
    if (value > maxV[ix]) maxV[ix] = value;
  }
  bool isConstant(size_t ix) const { return minV[ix] == maxV[ix]; }
  bool hasSpread(size_t ix) const { return maxV[ix] > minV[ix]; }
  std::vector<double> minV, maxV;
};

// The squared deviations of values below about 1e-162 in magnitude
// underflow to 0. The column then keeps its values, which are too small to
// change a score.
void warnVarianceUnderflow(size_t featureIndex, const ColumnRanges& ranges,
                           size_t ix) {
  const std::streamsize precision = std::cerr.precision(6);
  std::cerr << "Warning: the variance of feature "
            << DataSet::getFeatureNames().describeFeature(featureIndex)
            << " underflows to 0, as its values, between " << ranges.minV[ix]
            << " and " << ranges.maxV[ix] << ", are too small in magnitude. "
            << "The feature is not scaled and has practically no effect on "
            << "the scores. Rescale it in the input, e.g. multiply it by a "
            << "constant." << std::endl;
  std::cerr.precision(precision);
}

}  // namespace

StdvNormalizer::StdvNormalizer() {
}

StdvNormalizer::~StdvNormalizer() {
}

void StdvNormalizer::unnormalizeweight(const std::vector<double>& in, 
    std::vector<double>& out) const {
  double sum = 0;
  unsigned int i = 0;
  for (; i < numFeatures; i++) {
    out[i] = in[i] / div[i];
    sum += sub[i] * in[i] / div[i];
  }
  out[i] = in[i] - sum;
}

void StdvNormalizer::normalizeweight(const std::vector<double>& in, 
    std::vector<double>& out) const {
  double sum = 0;
  size_t i = 0;
  for (; i < numFeatures; i++) {
    out[i] = in[i] * div[i];
    sum += sub[i] * in[i];
  }
  out[i] = in[i] + sum;
}

void StdvNormalizer::setSet(std::vector<double*>& featuresV,
                            std::vector<double*>& rtFeaturesV, size_t nf,
                            size_t nrf) {
  numFeatures = nf;
  numRetentionFeatures = nrf;
  sub.resize(nf + nrf, 0.0);
  div.resize(nf + nrf, 0.0);
  double n = 0.0;
  double* features;
  size_t ix;
  vector<double*>::iterator it = featuresV.begin();
  ColumnRanges ranges(nf + nrf);
  for (; it != featuresV.end(); ++it) {
    features = *it;
    n++;
    for (ix = 0; ix < numFeatures; ++ix) {
      sub[ix] += features[ix];
      ranges.add(ix, features[ix]);
    }
  }
  for (it = rtFeaturesV.begin(); it != rtFeaturesV.end(); ++it) {
    features = *it;
    for (ix = numFeatures; ix < numFeatures + numRetentionFeatures; ++ix) {
      sub[ix] += features[ix - numFeatures];
      ranges.add(ix, features[ix - numFeatures]);
    }
  }
  // the factors are printed with two digits; the precision of cerr is put
  // back afterwards, since the Stage lines of -v 3 are printed with it
  const std::streamsize oldPrecision = cerr.precision();
  if (VERB > 2) {
    cerr.precision(2);
    cerr << "Normalization factors" << endl << "Avg ";
  }
  for (ix = 0; ix < numFeatures + numRetentionFeatures; ++ix) {
    if (n > 0.0) {
      sub[ix] /= n;
    }
    // The sum of n copies of a value such as 0.1 is rounded, so the mean of
    // a constant column can differ from its value in the last bits. Every
    // row then has the same tiny deviation, which the standard
    // deviation scales to 1 or -1: a second bias column instead of a
    // column of zeros. The mean of a constant column is its value.
    if (ranges.isConstant(ix)) {
      sub[ix] = ranges.minV[ix];
    }
    if (VERB > 2) {
      cerr << "\t" << sub[ix];
    }
  }
  for (it = featuresV.begin(); it != featuresV.end(); ++it) {
    features = *it;
    for (ix = 0; ix < numFeatures; ++ix) {
      if (!std::isfinite(features[ix])) {
        cerr << "Reached strange feature with val=" << features[ix]
            << " at col=" << ix << endl;
      }
      double d = features[ix] - sub[ix];
      div[ix] += d * d;
    }
  }
  for (it = rtFeaturesV.begin(); it != rtFeaturesV.end(); ++it) {
    features = *it;
    for (ix = numFeatures; ix < numFeatures + numRetentionFeatures; ++ix) {
      if (!std::isfinite(features[ix-numFeatures])) {
        cerr << "Reached strange feature with val=" << features[ix
            - numFeatures] << " at col=" << ix << endl;
      }
      double d = features[ix - numFeatures] - sub[ix];
      div[ix] += d * d;
    }
  }
  if (VERB > 2) {
    cerr << endl << "Stdv";
  }
  for (ix = 0; ix < numFeatures + numRetentionFeatures; ++ix) {
    if (div[ix] <= 0 || n == 0) {
      if (ranges.hasSpread(ix)) {
        warnVarianceUnderflow(ix, ranges, ix);
      }
      div[ix] = 1.0;
    } else {
      div[ix] = sqrt(div[ix] / n);
    }
    if (VERB > 2) {
      cerr << "\t" << div[ix];
    }
  }
  if (VERB > 2) {
    cerr << endl;
    cerr.precision(oldPrecision);
  }
  checkFactorsFinite(0, numFeatures + numRetentionFeatures, "mean",
                     "standard deviation");
}

}}}  // namespace OpenMS::Internal::Percolator
