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

#include "FeatureMemoryPool.h"

void FeatureMemoryPool::createPool(size_t numFeatures) {
  numFeatures_ = static_cast<unsigned int>(numFeatures);
  numRowsPerBlock_ = kBlockSize / numFeatures_;
  isInitialized_ = true;
}

void FeatureMemoryPool::createNewBlock() {
  double* memStart = new double[numFeatures_ * numRowsPerBlock_]();
  memStarts_.push_back(memStart);
}

void FeatureMemoryPool::destroyPool() {
  for (size_t i = 0; i < memStarts_.size(); ++i) {
    if (memStarts_.at(i) != NULL) {
      delete[] memStarts_.at(i);
      memStarts_.at(i) = NULL;
    }
  }
  isInitialized_ = false;
}

double* FeatureMemoryPool::addressFromIdx(unsigned int i) const {
  return memStarts_.at(i / numRowsPerBlock_) + (i % numRowsPerBlock_) * numFeatures_;
}

double* FeatureMemoryPool::allocate() {
  if (freeRows_.size() == 0) {
    if (initializedRows_ >= numRowsPerBlock_ * memStarts_.size()) {
      createNewBlock();
    }
    freeRows_.push_back(addressFromIdx(initializedRows_));
    initializedRows_++;
  }
  double* ret = freeRows_.back();
  freeRows_.pop_back();
  return ret;
}

void FeatureMemoryPool::deallocate(double* p) {
  freeRows_.push_back(p);
}
