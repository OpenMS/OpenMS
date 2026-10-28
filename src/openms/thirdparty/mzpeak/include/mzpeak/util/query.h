/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <any>
#include <boost/compat/function_ref.hpp>
#include <functional>
#include <optional>
#include <utility>
#include <variant>

#include "mzpeak/schema/group.h"
#include "mzpeak/util/types.h"

namespace MzPeak::Util {

/**
 * A low-level interface for selecting which records to extract from a
 * Parquet file.
 *
 * Queries are built using the Builder class and one of the provided
 * predicate functions.
 *
 * More complex queries can be constructed using the logic operators
 * (`&&`, `||`, and `!`) using the `and_then`, `or_else`, and `negate`
 * functions.
 */
class Query final {
private:
  struct Predicate;
  enum class Op { EQ, GT, LT, GE, LE };

public:
  /**
   * This class is used to construct a Query object using two inputs:
   *
   * 1. A Parquet group and field to compare to (`Schema::Column`)
   *
   * 2. A predicate function with a comparison value.
   *
   * NOTE: The C++ type given to the predicate functions (e.g., `eq`,
   * `gt`) must be compatible with the field type.  Currently this
   * is rather strict.  For example, if the field type is a PSI Int32,
   * then the predicate value *must* be an `int32_t`.
   *
   * Unfortunately this can only be checked at run-time without making
   * the query interface very difficult to use.  Therefore, mismatches
   * are reported as run-time exceptions.
   */
  class Builder final {
  public:
    /// Constructor.
    Builder(Schema::Column destination)
        : dest_(destination)
    {
    }

    /**
     * Field must match `val` exactly.
     */
    template <supported_type T> Query eq(T val) const
    {
      return validate({dest_, Op::EQ, val});
    }

    /**
     * Field must be greater than `val`.
     */
    template <supported_type T> Query gt(T val) const
    {
      return validate({dest_, Op::GT, val});
    }

    /**
     * Field must be less than `val`.
     */
    template <supported_type T> Query lt(T val) const
    {
      return validate({dest_, Op::LT, val});
    }

    /**
     * Field must be greater than or equal to `val`.
     */
    template <supported_type T> Query ge(T val) const
    {
      return validate({dest_, Op::GE, val});
    }

    /**
     * Field must be less than or equal to `val`.
     */
    template <supported_type T> Query le(T val) const
    {
      return validate({dest_, Op::LE, val});
    }

    /**
     * Low-level function for build a query with an operator.
     */
    template <supported_type T> Query via(T val, Op op) const
    {
      return validate({dest_, op, val});
    }

  private:
    Query validate(Predicate&& p) const;
    Schema::Column dest_;
  };

public:
  /**
   * A class used to return values to the query engine, and also the
   * final result return from query evaluation.
   *
   * This class models three possible states:
   *
   * 1. Failure.  The query should be terminated.
   *
   * 2. Null.  Treated like SQL NULL values.  That is, not an error
   * but might cause queries to short circuit.  Useful to signal that
   * certain columns can't be read because they are NULL.
   *
   * 3. Contains a valid value.
   */
  template <typename T> class Result {
  public:
    /// Unrecoverable failure.
    static Result fail() { return Result(true, std::nullopt); }

    /// NULL.
    static Result skip() { return Result(false, std::nullopt); }

    /// Valid value.
    Result(T);

    /// Did the column request fail?
    bool failed() const { return failed_; }

    /// Not failed and not NULL.
    bool has_value() const;

    /// Get the actual value recorded.
    T value() const;

    /// Return true if the result is not failed, not skipped, has a
    /// value, and that value is the given value.
    bool is(T) const;

    /// Convert to another type while preserving failure and NULL
    /// status.  That is, U is ignored if the current result is NULL.
    template <typename U> Result<U> to(U) const;

    /// Combine values using logical operations.
    Result and_then(const Result& other);
    Result or_else(const Result& other);
    Result negate();

  private:
    explicit Result(bool f, std::optional<T> r)
        : failed_(f)
        , result_(r)
    {
    }

    // This is so stupid.
    template <typename U> friend class Result;

    bool failed_;
    std::optional<T> result_;
  };

public:
  /// Join two queries together with a logical AND.
  Query and_then(const Query&) const;

  /// Join two queries together with a logical OR.
  Query or_else(const Query&) const;

  /// Negate a query.
  Query negate() const;

  // Column types that can be used in a query.
  using value_t = any_value_type; // From types.h

  // Like value_t, but instead of a single value this type is used for
  // evaluating a query on a min/max range.
  using range_t = any_pair_type; // From types.h

  /// A function that when given an column type, should return a single value.
  /// If this isn't possible it should return `Result<value_t>::skip()`.
  using eval_callback_t =
      boost::compat::function_ref<Result<value_t>(Schema::Column)>;

  /// A func ion that when given an column type should return a min
  /// and max.  If this isn't possible it should return
  /// `Result<range_t>::skip()`.
  using eval_range_callback_t =
      boost::compat::function_ref<Result<range_t>(Schema::Column)>;

  /**
   * Evaluate a query.
   */
  Result<bool> eval(eval_callback_t) const;

  /**
   * When this query is exactly one equality test on one column, the column and
   * the value it must equal; otherwise nothing.
   *
   * Selecting one entity by index is overwhelmingly the query this library
   * runs, and evaluating a general predicate tree per row -- through a
   * std::function, a type dispatch and two shared_ptr copies each time -- costs
   * far more than the comparison it performs.  This lets the executor resolve
   * the column once and loop over raw values instead.
   */
  std::optional<std::pair<Schema::Column, value_t>> as_equality() const;

  /**
   * Evaluate a range query.
   *
   * Range queries test to see if the query would match a value within
   * a min or max range.
   */
  Result<bool> eval(eval_range_callback_t) const;

private:
  friend class Builder;
  template <typename Fn, typename V> friend struct EvalHelper;

  // A predicate that can be tested against a value.
  struct Predicate {
    Schema::Column dest;
    Op op;
    value_t val;
  };

  // A logical node with a connective.
  struct Node {
    enum class Connective { AND, OR };
    Connective connective_;
    std::any lhs_;
    std::any rhs_;
  };

  // A tree is a leaf or a node.
  using tree_t = std::variant<Predicate, Node>;

  // Private constructors.
  explicit Query(Predicate&& pred);
  explicit Query(Node&&);

  // Private helper function.
  Query join(const Query& other, Node::Connective conn) const;

  // Data members.
  tree_t tree_;
  bool not_ = false;
};

/******************************************************************************/
template <typename T>
Query::Result<T>::Result(T v)
    : failed_(false)
    , result_(v)
{
}

/******************************************************************************/
template <typename T> bool Query::Result<T>::has_value() const
{
  return !failed_ && result_.has_value();
}

/******************************************************************************/
template <typename T> T Query::Result<T>::value() const { return result_.value(); }

/******************************************************************************/
template <typename T> bool Query::Result<T>::is(T t) const
{
  if (failed_) return false;
  return has_value() && value() == t;
}

/******************************************************************************/
template <typename T>
template <typename U>
Query::Result<U> Query::Result<T>::to(U u) const
{
  std::optional<U> r;
  if (has_value()) r = u;
  return Result<U>(failed_, r);
}

/******************************************************************************/
template <typename T>
Query::Result<T> Query::Result<T>::and_then(const Query::Result<T>& other)
{
  if (failed_) return *this;
  if (other.failed_) return other;
  if (!result_.has_value()) return *this;
  if (!other.result_.has_value()) return other;
  return Result(result_.value() && other.result_.value());
}

/******************************************************************************/
template <typename T>
Query::Result<T> Query::Result<T>::or_else(const Query::Result<T>& other)
{
  if (failed_) return *this;
  if (other.failed_) return other;
  if (!result_.has_value()) return other;
  if (!other.result_.has_value()) return *this;
  return Result(result_.value() || other.result_.value());
}

/******************************************************************************/
template <typename T> Query::Result<T> Query::Result<T>::negate()
{
  Result r = *this;
  r.result_ = r.result_.and_then([](auto& v) -> std::optional<T> { return !v; });
  return r;
}

} // namespace MzPeak::Util
