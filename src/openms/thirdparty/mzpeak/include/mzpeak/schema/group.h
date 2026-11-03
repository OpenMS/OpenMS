/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <string>

#include "mzpeak/schema/cv.h"
#include "mzpeak/schema/file.h"
#include "mzpeak/util/numpress.h"
#include "mzpeak/util/types.h"

// Forward declarations.
namespace parquet::schema {
class Node;
class GroupNode;
} // namespace parquet::schema

namespace MzPeak::Schema {

/**
 * Internal representation of a parquet schema `Group` with fields.
 *
 * The goal is to flatten the parquet schema as much as possible given
 * what we know about the mzPeak schema.
 */
class Group final {
public:
  /// Type used to store column indexes.
  using index_type = int32_t;

  // Type-safe wrapper for field CV types.
  struct CVType : CV {
    CVType(std::string code, std::string accession)
        : CV(std::move(code), std::move(accession))
    {
    }
  };

  // Type-safe wrapper for field CV units.
  struct CVUnit : CV {
    CVUnit(std::string code, std::string accession)
        : CV(std::move(code), std::move(accession))
    {
    }
  };

  /**
   * A possibly non-scalar field.
   */
  class Field final {
  public:
    enum class Kind {
      /// Elements are scalars values.
      Scalar,

      /// Elements are a list of scalars.
      List,

      /// Metadata parameters,
      Params,

      /// Elements are a struct of several named children, read as a whole.
      /// The grid-encoding column `<array>_grid` is one.
      Struct,

      /// No clue.
      Unknown,
    };

    /// Constructor from an encoded column name.
    explicit Field(std::string_view column_name,
                   index_type rel_index,
                   index_type abs_index);

    /// Destructor.
    ~Field() = default;

    /**
     * Column index inside the parent group.
     */
    index_type relative_index() const;

    /**
     * Column index inside the parquet file.
     */
    index_type absolute_index() const;

    /**
     * The name of this field as recognized by parquet.
     */
    const std::string& name() const;

    /**
     * The structural type this field represents.
     */
    Kind kind() const;

    /**
     * Controlled vocabulary code and accession for the field type.
     */
    const std::optional<CVType>& cv_type() const;

    /**
     * Controlled vocabulary code and accession for the field unit.
     */
    const std::optional<CVUnit>& cv_unit() const;

    /**
     * The data type for values in this field.
     */
    const std::optional<Util::Type>& type() const;

    /**
     * Update the field's data type.
     */
    void type(Util::Type);

    /**
     * Return the numpress method type if the column name indicates
     * this is a numpress compressed column of `uint8_t`.
     */
    std::optional<Util::Numpress::Type> possibly_numpress() const;

  private:
    friend class Group;

    index_type rel_index_;
    index_type abs_index_;
    std::string schema_name_;
    std::optional<CVType> cv_type_;
    std::optional<CVUnit> cv_unit_;
    std::optional<Util::Type> type_;
    Kind kind_ = Kind::Scalar;
  };

public:
  /// Fields are stored in a map for quick look-up using their name.
  using field_map_t = std::map<std::string, std::shared_ptr<Field>>;

  /// Constructor for the root group to hold all of the top-level
  /// columns that are not in a separate struct/group.
  explicit Group(const parquet::schema::GroupNode&, const Schema::File&);

  /// Constructor from a parquet schema descriptor.
  explicit Group(const parquet::schema::GroupNode&,
                 const Schema::File&,
                 index_type index,
                 index_type offset);

  /// Destructor.
  ~Group() = default;

  /**
   * The schema name for this group.
   */
  const std::string& name() const;

  /**
   * Return the schema column index of this group.
   */
  index_type index() const;

  /**
   * Find a field given its name.
   *
   * NOTE: For metadata groups this is the cleaned name, not the raw
   * schema node name.
   */
  std::optional<std::shared_ptr<const Field>> field(std::string_view) const;

  /**
   * Find a field given its CV type.
   */
  std::optional<std::shared_ptr<const Field>> field(const CVType&&) const;

  /**
   * Return a map of all fields.
   */
  const field_map_t& fields() const;

  /**
   * Return true if this is the root group.
   *
   * There is only one root group and its fields represent the
   * top-level columns that are not themselves members of a group or
   * struct.
   */
  bool is_root() const;

  /**
   * Return a schema path to the given field.  Mostly useful for error
   * messages.
   */
  std::string path(const Field&) const;

private:
  void
  make_fields(const parquet::schema::GroupNode&, const Schema::File&, index_type);

  std::string name_;
  bool is_root_;
  index_type index_;
  field_map_t fields_;
};

/**
 * Mapping from group name to a Group.
 */
using GroupMap = std::map<std::string, std::shared_ptr<Group>>;

/**
 * A parquet column can be uniquely identified using its parent group
 * and field.
 */
using Column =
    std::pair<std::shared_ptr<const Group>, std::shared_ptr<const Group::Field>>;

} // namespace MzPeak::Schema
