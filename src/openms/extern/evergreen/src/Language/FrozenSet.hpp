#ifndef _FROZENSET_HPP
#define _FROZENSET_HPP

// From Python's frozenset class:
// See: https://stackoverflow.com/questions/20832279/python-frozenset-hashing-algorithm-implementation
template <typename K>
struct SetFrozenHash {
  std::size_t operator() (const std::set<K> & s) const {
    std::hash<K> single_hash;
    std::size_t combined_hash_value = 0;

    for (const K & obj : s) {
      unsigned long single = single_hash(obj);
      combined_hash_value ^= (single^(single<<16)^89869747ul) * 3644798167ul;
    }
    return combined_hash_value * 69069 + 907133923;
  }
};

template <typename K>
class FrozenSet {
private:
  std::set<K> _data;
  std::size_t _hash_value;

public:
  FrozenSet(const std::set<K> & s):
    _data(s)
  {
    SetFrozenHash<K> sh;
    _hash_value = sh(_data);
  }

  const std::set<K> & get_set() const {
    return _data;
  }
  const std::size_t hash_value() const {
    return _hash_value;
  }

  friend bool operator <(const FrozenSet<K> & lhs, const FrozenSet<K> & rhs) {
    return lhs._data < rhs._data;
  }

  friend bool operator ==(const FrozenSet<K> & lhs, const FrozenSet<K> & rhs) {
    return lhs._data == rhs._data;
  }
};

template <typename K>
  struct FrozenSetHash {
    std::size_t operator() (const FrozenSet<K> & fs) const {
      return fs.hash_value();
  }
};
#endif
