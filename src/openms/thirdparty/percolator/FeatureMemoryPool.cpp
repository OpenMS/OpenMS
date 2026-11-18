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

#include <cstdlib>
#include <new>

#if defined(__linux__)
#include <sys/mman.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {


namespace {
#if defined(__linux__)
const std::size_t kHugePageSize = std::size_t(1) << 21;
#endif
}  // namespace

void FeatureMemoryPool::BlockDeleter::operator()(double* block) const {
#if defined(__linux__)
  free(block);
#else
  delete[] block;
#endif
}

void FeatureMemoryPool::createPool(size_t numFeatures) {
  numFeatures_ = static_cast<unsigned int>(numFeatures);
  rowStride_ = numFeatures_ + (svmRowsInPlace_ ? 1u : 0u);
  numRowsPerBlock_ = kBlockSize / rowStride_;
  isInitialized_ = true;
}

void FeatureMemoryPool::createNewBlock() {
  // every row is completely written by its user, so no zero-initialization
  const std::size_t numDoubles =
      static_cast<std::size_t>(rowStride_) * numRowsPerBlock_;
#if defined(__linux__)
  void* block = NULL;
  if (hasInterceptColumn()) {
    // With in-place training the SVM reads the rows in score order, i.e.
    // scattered over the pool, and 4 KiB pages would cost a TLB miss for
    // most rows. Blocks aligned to 2 MiB can be backed by transparent huge
    // pages, also when the system enables them only on request (madvise).
    const std::size_t bytes =
        (numDoubles * sizeof(double) + kHugePageSize - 1) / kHugePageSize *
        kHugePageSize;
    if (posix_memalign(&block, kHugePageSize, bytes) != 0) {
      block = NULL;
    }
#ifdef MADV_HUGEPAGE
    if (block != NULL) {
      madvise(block, bytes, MADV_HUGEPAGE);  // a hint, failure is harmless
    }
#endif
  } else {
    block = malloc(numDoubles * sizeof(double));
  }
  if (block == NULL) {
    throw std::bad_alloc();
  }
  memStarts_.emplace_back(static_cast<double*>(block));
#else
  memStarts_.emplace_back(new double[numDoubles]);
#endif
}

void FeatureMemoryPool::destroyPool() {
  memStarts_.clear();
  freeRows_.clear();
  initializedRows_ = 0;
  isInitialized_ = false;
}

double* FeatureMemoryPool::addressFromIdx(unsigned int i) const {
  return memStarts_.at(i / numRowsPerBlock_).get() +
         static_cast<std::size_t>(i % numRowsPerBlock_) * rowStride_;
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
  if (hasInterceptColumn()) {
    ret[numFeatures_] = 1.0;  // the intercept column of the SVM
  }
  return ret;
}

void FeatureMemoryPool::deallocate(double* p) {
  freeRows_.push_back(p);
}

}}}  // namespace OpenMS::Internal::Percolator
