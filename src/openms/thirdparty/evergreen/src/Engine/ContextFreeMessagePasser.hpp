#ifndef _CONTEXTFREEMESSAGEPASSER_HPP
#define _CONTEXTFREEMESSAGEPASSER_HPP

#include "MessagePasser.hpp"

// ContextFreeMessagePasser types do not have distinct roles for
// messages in or out; this contrasts with e.g. probabilistic
// addition, where the output edge needs to be labeled as distinct
// from the input edges. As a result, ContextFreeMessagePasser types
// can add edges outside of their constructors; therefore, promote
// bind_to to allow public access.

// ContextFreeMessagePasser is a good base type for hyperedges.

template <typename VARIABLE_KEY>
class ContextFreeMessagePasser : public MessagePasser<VARIABLE_KEY> {
protected:
  ContextFreeMessagePasser()
  { }

public:
  using MessagePasser<VARIABLE_KEY>::bind_to;
};

#endif
