/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include "mzpeak/util/query.h"

#include <cassert>
#include <functional>
#include <variant>

#include "mzpeak/exception.h"
#include "mzpeak/util/compat.h" // IWYU pragma: keep

namespace MzPeak::Util {

// There used to be a Trampoline<T> here, "to turn recursive algorithms into
// iteration". It never did: no caller ever built a thunk, every Trampoline was
// constructed from an already-computed Result, and eval/eval_node recursed
// directly through it. Its one observable effect was a compile error on
// Apple's clang 15-17, whose libc++ lacks std::move_only_function and falls
// back to Boost.Compat's, which checks invocability of the thunk's
// `Trampoline<T>()` signature while Trampoline<T> is still being defined.
// Query trees are a few levels deep; plain recursion is the honest shape.

/******************************************************************************/
std::pair<Query::value_t, Query::value_t> decode_range_type(const Query::range_t& rt)
{
  return std::visit(
      [](auto&& pair) {
        auto first = Query::value_t{pair.first};
        auto second = Query::value_t{pair.second};
        return std::make_pair(first, second);
      },
      rt);
}

/******************************************************************************/
Query Query::Builder::validate(Query::Predicate&& p) const
{
  if (!p.dest.second->type().has_value()) {
    std::string msg("cannot query field with unknown type: ");
    msg += p.dest.first->path(*p.dest.second);
    throw TypeError(msg);
  }

  std::visit([&type = p.dest.second->type().value()](
                 auto&& v) -> void { ensure_same_type<decltype(v)>(type); },
             p.val);

  return Query(std::move(p));
}

/******************************************************************************/
template <typename Fn, typename V> struct EvalHelper {
  using result_t = Query::Result<bool>;

  // Eval a query using the given function for fetching values.
  Query::Result<bool> eval(const Query& query, Fn fn) const;

  // Eval a query node.
  Query::Result<bool> eval_node(const Query::Node& node, Fn fn) const;

  // Dispatch on the type of the predicate's value.
  template <Util::Type T>
  Query::Result<bool> eval_predicate(const Query::Predicate& p, Fn fn) const;

  // Match a predicate against a single value.
  bool match(const Query::Predicate&, Query::value_t) const;

  // Match a predicate against a min/max range.
  bool match(const Query::Predicate& p, Query::range_t) const;
};

/******************************************************************************/
Query::Query(Predicate&& p)
    : tree_({std::move(p)})
{
}

/******************************************************************************/
Query::Query(Node&& c)
    : tree_({std::move(c)})
{
}

/******************************************************************************/
Query Query::negate() const
{
  Query q(*this);
  q.not_ = !q.not_;
  return q;
}

/******************************************************************************/
Query Query::and_then(const Query& rhs) const
{
  return join(rhs, Node::Connective::AND);
}

/******************************************************************************/
Query Query::or_else(const Query& rhs) const
{
  return join(rhs, Node::Connective::OR);
}

/******************************************************************************/
Query Query::join(const Query& other, Node::Connective oper) const
{
  return Query(Node{oper, *this, other});
}

/******************************************************************************/
Query::Result<bool> Query::eval(eval_callback_t fn) const
{
  EvalHelper<eval_callback_t, value_t> eh;
  return eh.eval(*this, fn);
}

/******************************************************************************/
Query::Result<bool> Query::eval(eval_range_callback_t fn) const
{
  EvalHelper<eval_range_callback_t, range_t> eh;
  return eh.eval(*this, fn);
}

/******************************************************************************/
template <typename Fn, typename V>
bool EvalHelper<Fn, V>::match(const Query::Predicate& p, Query::value_t v) const
{
  assert(v.index() == p.val.index());

  switch (p.op) {
  case Query::Op::EQ:
    return v == p.val;

  case Query::Op::GT:
    return v > p.val;

  case Query::Op::LT:
    return v < p.val;

  case Query::Op::GE:
    return v >= p.val;

  case Query::Op::LE:
    return v <= p.val;
  }

  return false;
}

/******************************************************************************/
template <typename Fn, typename V>
bool EvalHelper<Fn, V>::match(const Query::Predicate& p, Query::range_t v) const
{
  auto [min, max] = decode_range_type(v);
  assert(min.index() == max.index() && min.index() == p.val.index());

  switch (p.op) {
  case Query::Op::EQ:
    return (min == p.val || max == p.val) || (p.val > min && p.val < max);

  case Query::Op::GT:
    return max > p.val;

  case Query::Op::LT:
    return min < p.val;

  case Query::Op::GE:
    return max >= p.val;

  case Query::Op::LE:
    return min <= p.val;
  }

  return false;
}

/******************************************************************************/
template <typename Fn, typename V>
template <Util::Type T>
Query::Result<bool> EvalHelper<Fn, V>::eval_predicate(const Query::Predicate& p,
                                                      Fn fn) const
{
  Query::Result<V> val(std::invoke(fn, p.dest));
  if (!val.has_value()) return val.template to<bool>(false);

  return std::visit(
      [&](auto&& v) -> Query::Result<bool> {
        using T1 = Util::type_traits<T>::value_type;
        using T2 = std::decay_t<decltype(v)>;
        using T3 = std::pair<T1, T1>;

        if constexpr (std::is_same_v<T1, T2>) {
          return match(p, v);
        } else if constexpr (std::is_same_v<T2, T3>) {
          return match(p, v);
        } else {
          std::string msg("predicate and value mismatch: ");
          throw TypeError(msg);
          return Query::Result<bool>::fail(); // clang is too stupid to see the throw
        }
      },
      val.value());
}

/******************************************************************************/
template <typename Fn, typename V>
Query::Result<bool> EvalHelper<Fn, V>::eval(const Query& query, Fn fn) const
{
  result_t result = std::visit(
      [&](auto& tree) -> result_t {
        using T = std::decay_t<decltype(tree)>;

        if constexpr (std::is_same_v<T, Query::Predicate>) {
          if (!tree.dest.second->type().has_value()) {
            std::string msg("invalid query on field with unknown data type: ");
            msg += tree.dest.first->path(*tree.dest.second);
            throw TypeError(msg);
          }

          return lift_type(tree.dest.second->type().value(), [&]<Util::Type T> {
            return eval_predicate<T>(tree, fn);
          });
        } else if constexpr (std::is_same_v<T, Query::Node>) {
          return eval_node(tree, fn);
        } else {
          static_assert(false_type<T>, "unhanded variant");
        }
      },
      query.tree_);

  if (query.not_) {
    return result.negate();
  } else {
    return result;
  }
}

/******************************************************************************/
template <typename Fn, typename V>
Query::Result<bool> EvalHelper<Fn, V>::eval_node(const Query::Node& node, Fn fn) const
{
  result_t result = true;

  if (node.lhs_.has_value()) {
    result = eval(std::any_cast<Query>(node.lhs_), fn);
  }

  if (node.rhs_.has_value()) {
    switch (node.connective_) {
    case Query::Node::Connective::AND:
      if (result.is(true)) {
        return eval(std::any_cast<Query>(node.rhs_), fn);
      }
      break;
    case Query::Node::Connective::OR:
      if (result.is(false)) {
        return eval(std::any_cast<Query>(node.rhs_), fn);
      }
      break;
    }
  }

  return result;
}

/******************************************************************************/
std::optional<std::pair<Schema::Column, Query::value_t>> Query::as_equality() const
{
  // A negated equality is an inequality, and an inequality selects almost every
  // row -- the fast path below would be a pessimisation as well as wrong.
  if (not_) return std::nullopt;
  if (!std::holds_alternative<Predicate>(tree_)) return std::nullopt;

  const Predicate& predicate = std::get<Predicate>(tree_);
  if (predicate.op != Op::EQ) return std::nullopt;

  return std::make_pair(predicate.dest, predicate.val);
}

} // namespace MzPeak::Util
