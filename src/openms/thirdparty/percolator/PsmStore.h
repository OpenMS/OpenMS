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
#ifndef PSM_STORE_H_
#define PSM_STORE_H_

#include <cstring>
#include <deque>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "PSMDescription.h"

namespace OpenMS { namespace Internal { namespace Percolator {


/*
 * PsmStore owns PSMDescription objects in fixed-size blocks, together with
 * the records of their ids, peptides and protein lists. All addresses stay
 * valid while more PSMs are added, and there is no heap allocation per PSM.
 * Protein names are stored once each.
 */
class PsmStore {
 public:
  PsmStore() = default;
  PsmStore(const PsmStore&) = delete;
  PsmStore& operator=(const PsmStore&) = delete;

  // moves psm into the store and returns its address in the store
  PSMDescription* add(PSMDescription&& psm) {
    PSMDescription* slot;
    if (!freeSlots_.empty()) {
      slot = freeSlots_.back();
      freeSlots_.pop_back();
    } else {
      if (numUsed_ == blocks_.size() * kBlockSize) {
        blocks_.emplace_back(new PSMDescription[kBlockSize]);
      }
      slot = &blocks_.back()[numUsed_ % kBlockSize];
      ++numUsed_;
    }
    *slot = std::move(psm);
    return slot;
  }

  // empties the slot of a PSM so that a later add can reuse it
  void release(PSMDescription* psm) {
    *psm = PSMDescription();
    freeSlots_.push_back(psm);
  }

  // Copies the id, the peptide and the protein names of psm into the store,
  // each protein name only once, and points psm to the copies; see
  // PSMDescription::setStrings.
  void storeStrings(PSMDescription& psm, std::string_view id,
                    std::string_view peptide,
                    const std::vector<std::string_view>& proteins) {
    proteinList_.resize(proteins.size());
    for (size_t i = 0; i < proteins.size(); ++i) {
      auto found = proteinIndex_.find(proteins[i]);
      if (found == proteinIndex_.end()) {
        proteinNames_.emplace_back(proteins[i]);  // deque elements never move
        found = proteinIndex_
                    .emplace(std::string_view(proteinNames_.back()),
                             &proteinNames_.back())
                    .first;
      }
      proteinList_[i] = found->second;
    }
    const size_t size =
        PSMDescription::recordSize(id, peptide, proteins.size());
    psm.setStrings(allocate(size), id, peptide, proteinList_.data(),
                   proteinList_.size());
  }

  void clear() {
    blocks_.clear();
    freeSlots_.clear();
    numUsed_ = 0;
    charBlocks_.clear();
    longRecords_.clear();
    charUsed_ = 0;
    proteinIndex_.clear();
    proteinNames_.clear();
    proteinList_.clear();
  }

 private:
  static const size_t kBlockSize = 1u << 16;
  std::vector<std::unique_ptr<PSMDescription[]> > blocks_;
  std::vector<PSMDescription*> freeSlots_;
  size_t numUsed_ = 0;

  // records are packed into the last block; long ones get their own
  char* allocate(size_t size) {
    if (size > kCharBlockSize / 4) {  // too long to share a block
      longRecords_.emplace_back(new char[size]);
      return longRecords_.back().get();
    }
    if (charBlocks_.empty() || charUsed_ + size > kCharBlockSize) {
      charBlocks_.emplace_back(new char[kCharBlockSize]);
      charUsed_ = 0;
    }
    char* dest = charBlocks_.back().get() + charUsed_;
    charUsed_ += size;
    return dest;
  }

  static const size_t kCharBlockSize = 1u << 20;
  std::vector<std::unique_ptr<char[]> > charBlocks_, longRecords_;
  size_t charUsed_ = 0;
  std::deque<std::string> proteinNames_;
  std::unordered_map<std::string_view, const std::string*> proteinIndex_;
  std::vector<const std::string*> proteinList_;  // reused by storeStrings
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /* PSM_STORE_H_ */
