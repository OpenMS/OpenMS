#ifndef _SETHASH_HPP
#define _SETHASH_HPP

#include <functional>
#include <set>
#include <unordered_set>

template <typename T>
struct SetHash {
  // Mutlplying by a large prime should broadcast to higher bits (and
  // since it's prime, taking % 2^64 should be distributed fairly
  // uniformly). Also XOR with single element hash because product
  // with large prime may not use lower-significance bits effectively.
  std::size_t operator() (const std::unordered_set<T> & s) const {
    std::hash<T> single_hash;
    std::size_t combined_hash_value = 0;
    for (const T & obj : s) {
      unsigned long single = single_hash(obj);
      combined_hash_value += 2147483647ul*single ^ single;
    }
    combined_hash_value += 2147483647ul*s.size() ^ s.size();
    return combined_hash_value;
  }
};

#endif
