// Written by Oliver Serang 2009
// see license for more information

#ifndef _Set_H
#define _Set_H

#include "Array.h"

class Set : public Array<int> {
 public:
  Set() {}

  void add(int el);
  const Set & operator &=(const Set & rhs); // intersection
  const Set & operator |=(const Set & rhs); // union
  bool isEmpty() const { return size() == 0; }

  int find(int x) const;
  int findHelper(int low, int high, int x) const;

  Set reindexTo(const Set & rhs);
  Set reindexToFind(const Set & rhs);

  class UnsortedOrderException {};
  class InvalidBaseException {};

  Array<int> operator [](const Set & rhs) const {
    return Array<int>::operator [](rhs);
  }
  using Array<int>::operator[];
  const int & operator [](std::size_t k) const {
    return Array<int>::operator [](k);
  }

  // for hashing sets
  static unsigned int sumSetElements(const Set & s) {
    unsigned int sum = 0;
    for (std::size_t k=0; k<s.size(); k++)
      sum += static_cast<unsigned int>(s[k]);
    return sum;
  }

  using Array<int>::size;
  using Array<int>::Iterator;
  using Array<int>::begin;
  using Array<int>::end;

  Set without( const Set & rhs ) const;

  static Set FullSet(int low, int high);
  static Set SingletonSet(int value) {
    return FullSet(value, value);
  }

 private:
  bool verify() const;
};

Set operator &(Set lhs, const Set & rhs); // intersection
Set operator |(Set lhs, const Set & rhs); // union


#endif

