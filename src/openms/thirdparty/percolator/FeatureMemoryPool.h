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
#ifndef FEATURE_MEMORY_POOL_H_
#define FEATURE_MEMORY_POOL_H_

/* Adapted from https://www.thinkmind.org/download.php?articleid=computation_tools_2012_1_10_80006 */

#include <memory>
#include <vector>
#include <iostream>

namespace OpenMS { namespace Internal { namespace Percolator {


/* Storage of the feature rows of the PSMs, in blocks of 8 MiB. */
class FeatureMemoryPool {
 private:
   static const unsigned int kBlockSize = 1u << 20; // in number of doubles, e.g. 8MB if sizeof(double) = 8
   struct BlockDeleter {
     void operator()(double* block) const;
   };
   unsigned int numRowsPerBlock_, numFeatures_, rowStride_, initializedRows_;
   std::vector<std::unique_ptr<double[], BlockDeleter> > memStarts_;
   std::vector<double*> freeRows_;
   bool isInitialized_;
   bool svmRowsInPlace_;
 public:
  FeatureMemoryPool() : numRowsPerBlock_(0), numFeatures_(0), rowStride_(0),
                        initializedRows_(0), isInitialized_(false),
                        svmRowsInPlace_(false) {}

  FeatureMemoryPool(const FeatureMemoryPool&) = delete;
  FeatureMemoryPool& operator=(const FeatureMemoryPool&) = delete;

  void createPool(size_t numFeatures);
  void createNewBlock();
  void destroyPool();
  
  inline bool isInitialized() const { return isInitialized_; }

  // With on, the pools created afterwards give rows of numFeatures features
  // followed by the constant 1.0, the intercept column of the SVM, and put
  // their blocks on huge pages where the system allows it (Linux), so that
  // the SVM can train on the rows in place (see AlgIn and
  // CrossValidation::setSvmRowsInPlace).
  inline void setSvmRowsInPlace(bool on) { svmRowsInPlace_ = on; }
  inline bool svmRowsInPlace() const { return svmRowsInPlace_; }
  // number of doubles from one row to the next
  inline unsigned int rowStride() const { return rowStride_; }
  // whether the rows of the current pool end with 1.0
  inline bool hasInterceptColumn() const { return rowStride_ > numFeatures_; }
  // whether the rows of the current pool hold numFeatures features followed
  // by the intercept 1.0, as training in place needs
  inline bool canTrainInPlace(size_t numFeatures) const {
    return hasInterceptColumn() && numFeatures_ == numFeatures;
  }

  double* addressFromIdx(unsigned int i) const;

  double* allocate();
  void deallocate(double* p);
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /* FEATURE_MEMORY_POOL_H_ */
