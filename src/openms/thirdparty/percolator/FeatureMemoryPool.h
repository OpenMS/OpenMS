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

#include <vector>
#include <iostream>

class FeatureMemoryPool {
 private:
   static const unsigned int kBlockSize = 65536; // in number of doubles, e.g. 0.5MB if sizeof(double) = 8
   unsigned int numRowsPerBlock_, numFeatures_, initializedRows_;
   std::vector<double*> memStarts_;
   std::vector<double*> freeRows_;
   bool isInitialized_;
 public:
  FeatureMemoryPool() : numRowsPerBlock_(0), numFeatures_(0), 
                        initializedRows_(0), isInitialized_(false) {}

  ~FeatureMemoryPool() { destroyPool(); }

  void createPool(size_t numFeatures);
  void createNewBlock();
  void destroyPool();
  
  inline bool isInitialized() const { return isInitialized_; }

  double* addressFromIdx(unsigned int i) const;

  double* allocate();
  void deallocate(double* p);
};

#endif /* FEATURE_MEMORY_POOL_H_ */
