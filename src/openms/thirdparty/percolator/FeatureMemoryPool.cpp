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

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

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

/**
 * OpenMS: creates the blocks of the first numRows rows (if the pool has fewer), zero-filled as
 * createNewBlock does, by the threads of an OpenMP parallel loop: allocate() then hands out the
 * same rows in the same order without creating blocks (filling a block of a large set on one
 * thread took most of the time of allocating its rows).
 */
void FeatureMemoryPool::reserveRows(std::size_t numRows) {
  if (numRowsPerBlock_ == 0) return;
  const std::size_t numBlocks = (numRows + numRowsPerBlock_ - 1) / numRowsPerBlock_;
  if (numBlocks <= memStarts_.size()) return;
  const std::size_t blockSize = static_cast<std::size_t>(numFeatures_) * numRowsPerBlock_;
  std::vector<std::unique_ptr<double[]> > blocks(numBlocks - memStarts_.size());
  for (std::unique_ptr<double[]>& block : blocks) block.reset(new double[blockSize]);
  memStarts_.reserve(numBlocks);
#pragma omp parallel for schedule(static)
  for (std::ptrdiff_t b = 0; b < static_cast<std::ptrdiff_t>(blocks.size()); ++b) {
    std::fill_n(blocks[static_cast<std::size_t>(b)].get(), blockSize, 0.0);
  }
  for (std::unique_ptr<double[]>& block : blocks) memStarts_.push_back(block.release());
}

// OpenMS: the smallest block of rows a task of moveRowsInOrder copies
[[maybe_unused]] constexpr std::ptrdiff_t kMinRowsPerTask = 1024;

// OpenMS: whether moveRowsInOrder runs in an OpenMP parallel region (its copies as tasks)
static bool inParallel() {
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

/**
 * OpenMS: if *rows[0], *rows[1], ... are the rows of the pool (its initialized rows, none of
 * them free), each once, and rows have numFeatures elements (as in the pool): moves the row
 * *rows[k] to row k (addressFromIdx(k)), sets *rows[k] to it and returns true. Otherwise
 * returns false and changes nothing. (Moving the rows one
 * after the other by swaps, as Scores::createXvalSetsBySpectrum does otherwise, leaves the
 * same rows and pointers.)
 */
bool FeatureMemoryPool::moveRowsInOrder(const std::vector<double**>& rows,
                                        std::size_t numFeatures) {
  const std::size_t n = rows.size();
  if (!freeRows_.empty() || n != initializedRows_ || numFeatures_ == 0 ||
      numFeatures != numFeatures_) {
    return false;
  }
  // the blocks by address, to find the index of a row
  std::vector<std::pair<std::uintptr_t, std::size_t> > blocks;
  for (std::size_t b = 0; b < memStarts_.size(); ++b) {
    blocks.emplace_back(reinterpret_cast<std::uintptr_t>(memStarts_[b]), b);
  }
  std::sort(blocks.begin(), blocks.end());
  const std::uintptr_t rowBytes = static_cast<std::uintptr_t>(numFeatures_) * sizeof(double);
  std::vector<char> seen(n, 0);
  for (std::size_t k = 0; k < n; ++k) {
    const std::uintptr_t a = reinterpret_cast<std::uintptr_t>(*rows[k]);
    auto it = std::upper_bound(blocks.begin(), blocks.end(),
                               std::make_pair(a, std::numeric_limits<std::size_t>::max()));
    if (it == blocks.begin()) return false;
    --it;
    const std::uintptr_t offset = a - it->first;
    if (offset % rowBytes != 0 || offset / rowBytes >= numRowsPerBlock_) return false;
    const std::size_t idx = it->second * numRowsPerBlock_ + static_cast<std::size_t>(offset / rowBytes);
    if (idx >= n || seen[idx]) return false;
    seen[idx] = 1;
  }
  // copied out and back in blocks of rows, by the tasks of an enclosing OpenMP parallel region
  // if there is one (the blocks are disjoint)
  const std::unique_ptr<double[]> moved(new double[n * numFeatures_]);
  const std::ptrdiff_t numRows = static_cast<std::ptrdiff_t>(n);
  const bool inTasks = n >= 2 * kMinRowsPerTask && inParallel();
  (void)inTasks;
#pragma omp taskloop grainsize(kMinRowsPerTask) shared(rows, moved) if(inTasks)
  for (std::ptrdiff_t k = 0; k < numRows; ++k) {
    std::copy(*rows[k], *rows[k] + numFeatures_, moved.get() + k * numFeatures_);
  }
#pragma omp taskloop grainsize(kMinRowsPerTask) shared(rows, moved) if(inTasks)
  for (std::ptrdiff_t k = 0; k < numRows; ++k) {
    double* row = addressFromIdx(static_cast<unsigned int>(k));
    std::copy(moved.get() + k * numFeatures_, moved.get() + (k + 1) * numFeatures_, row);
    *rows[k] = row;
  }
  return true;
}

}}}  // namespace OpenMS::Internal::Percolator
