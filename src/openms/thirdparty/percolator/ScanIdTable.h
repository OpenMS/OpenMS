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
#ifndef SCAN_ID_TABLE_H_
#define SCAN_ID_TABLE_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace OpenMS { namespace Internal { namespace Percolator {


/*
 * Tells whether a scan id, i.e. a scan number and an experimental mass,
 * occurs both in a PSM with a decoy label and in a PSM with another label,
 * which means that the targets and decoys come from separate searches.
 *
 * Scan ids are equal when both the scan numbers and the masses compare
 * equal, as in a hash table keyed on the pair: 0.0 and -0.0 are the same
 * mass, and a NaN mass is never equal to another one.
 *
 * The scan ids are collected in a vector, which costs 16 bytes per PSM
 * instead of a hash table node, and sorted. The vector is sorted and checked
 * whenever its size doubles, by sorting the new part and merging it into the
 * sorted part. Separate searches therefore usually show up early, after which
 * the vector is released and further scan ids are ignored.
 */
class SeparateSearchDetector {
 public:
  // returns true once a scan id with both labels has been found
  bool add(int scan, double mass, bool isDecoy) {
    if (found_) return true;
    if (mass != mass) return false;  // a NaN mass matches no other scan id
    keys_.push_back(Key{scan, isDecoy, mass});
    if (keys_.size() >= nextCheck_) {
      check();
      nextCheck_ = 2 * keys_.size();
    }
    return found_;
  }

  // checks the scan ids added since the last check and releases the memory
  bool foundSeparateSearches() {
    if (!found_ && sorted_ < keys_.size()) check();
    std::vector<Key>().swap(keys_);
    sorted_ = 0;
    return found_;
  }

 private:
  struct Key {
    int scan;
    bool isDecoy;
    double mass;
  };
  // Orders by scan id, then by label. A total order on scan ids without
  // NaN masses: 0.0 and -0.0 are equivalent, as they are equal.
  static bool less(const Key& a, const Key& b) {
    if (a.scan != b.scan) return a.scan < b.scan;
    if (a.mass != b.mass) return a.mass < b.mass;
    return a.isDecoy < b.isDecoy;
  }

  void check() {
    std::vector<Key>::iterator mid = keys_.begin() + static_cast<std::ptrdiff_t>(sorted_);
    std::sort(mid, keys_.end(), less);
    std::inplace_merge(keys_.begin(), mid, keys_.end(), less);
    sorted_ = keys_.size();
    // within a run of equal scan ids, a decoy follows a non-decoy
    for (size_t i = 1; i < keys_.size(); ++i) {
      if (keys_[i].isDecoy != keys_[i - 1].isDecoy &&
          keys_[i].scan == keys_[i - 1].scan &&
          keys_[i].mass == keys_[i - 1].mass) {
        found_ = true;
        std::vector<Key>().swap(keys_);
        sorted_ = 0;
        return;
      }
    }
  }

  std::vector<Key> keys_;
  size_t sorted_ = 0;  // keys_[0, sorted_) is sorted
  size_t nextCheck_ = 1u << 16;
  bool found_ = false;
};

/*
 * Hash table from a scan id to the sampling priority of the first PSM with
 * that scan id, for the reservoir sampling of -N, and one byte that the
 * caller may use, e.g. for the label of that PSM. Scan ids are equal as in
 * SeparateSearchDetector.
 *
 * Open addressing with linear probing in one array of 16-byte entries plus
 * the bytes in a second array: no allocation per scan id, and less than half
 * the memory of a node-based hash table.
 */
class ScanIdPriorityTable {
 public:
  // If the scan id was stored before, sets found to true and priority to its
  // priority. Otherwise stores it with the priority newPriority() and the
  // byte 0, and sets found to false. Returns the byte of the scan id, which
  // stays valid until the next call. A scan id with a NaN mass is never
  // stored, so it gets a new priority and a new byte every time.
  template <class NewPriority>
  uint8_t& findOrInsert(int scan, double mass, NewPriority newPriority,
                        uint32_t& priority, bool& found) {
    if (mass != mass) {
      priority = static_cast<uint32_t>(newPriority());
      found = false;
      unstoredByte_ = 0u;
      return unstoredByte_;
    }
    if (10 * (size_ + 1) > 7 * slots_.size()) grow();
    size_t i = slotOf(scan, mass);
    for (; !isEmpty(slots_[i]); i = (i + 1) & mask_) {
      if (slots_[i].scan == scan && slots_[i].mass == mass) {
        priority = slots_[i].priority;
        found = true;
        return bytes_[i];
      }
    }
    priority = static_cast<uint32_t>(newPriority());
    slots_[i] = Slot{mass, scan, priority};
    bytes_[i] = 0u;
    ++size_;
    found = false;
    return bytes_[i];
  }

  size_t size() const { return size_; }

 private:
  struct Slot {
    double mass;  // NaN in an empty slot
    int scan;
    uint32_t priority;
  };
  static bool isEmpty(const Slot& slot) { return slot.mass != slot.mass; }

  size_t slotOf(int scan, double mass) const {
    if (mass == 0.0) mass = 0.0;  // the same hash for -0.0
    uint64_t bits;
    std::memcpy(&bits, &mass, sizeof(bits));
    uint64_t h = bits ^ (static_cast<uint64_t>(static_cast<uint32_t>(scan)) *
                         0x9E3779B97F4A7C15ULL);
    // the finalizer of MurmurHash3
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return static_cast<size_t>(h) & mask_;
  }

  void grow() {
    const Slot empty = {std::numeric_limits<double>::quiet_NaN(), 0, 0u};
    std::vector<Slot> oldSlots(std::max<size_t>(slots_.size() * 2, 1u << 16),
                               empty);
    std::vector<uint8_t> oldBytes(oldSlots.size(), 0u);
    oldSlots.swap(slots_);
    oldBytes.swap(bytes_);
    mask_ = slots_.size() - 1;
    for (size_t j = 0; j < oldSlots.size(); ++j) {
      if (isEmpty(oldSlots[j])) continue;
      size_t i = slotOf(oldSlots[j].scan, oldSlots[j].mass);
      while (!isEmpty(slots_[i])) i = (i + 1) & mask_;
      slots_[i] = oldSlots[j];
      bytes_[i] = oldBytes[j];
    }
  }

  std::vector<Slot> slots_;  // the size is a power of 2
  std::vector<uint8_t> bytes_;  // the byte of the scan id in slots_[i]
  size_t size_ = 0;
  size_t mask_ = 0;
  uint8_t unstoredByte_ = 0u;  // the byte of a scan id with a NaN mass
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /* SCAN_ID_TABLE_H_ */
