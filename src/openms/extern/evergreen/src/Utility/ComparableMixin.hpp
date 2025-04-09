#ifndef _COMPARABLEMIXIN_H
#define _COMPARABLEMIXIN_H

template <typename T>
class ComparableMixin {
public:
  friend bool operator >(T lhs, T rhs) {
    return !( lhs < rhs || lhs == rhs );
  }
  friend bool operator !=(T lhs, T rhs) {
    return !( lhs == rhs );
  }
  friend bool operator >=(T lhs, T rhs) {
    return lhs > rhs || lhs == rhs;
  }
  friend bool operator <=(T lhs, T rhs) {
    return lhs < rhs || lhs == rhs;
  }
};

#endif

