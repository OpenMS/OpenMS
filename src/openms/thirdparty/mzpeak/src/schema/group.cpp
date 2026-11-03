/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include <memory>
#include <parquet/schema.h>
#include <parquet/types.h>
#include <ranges>

#include "mzpeak/schema/group.h"
#include "mzpeak/schema/psi/data_type.h"

namespace MzPeak::Schema {

/******************************************************************************/
template <typename T> std::optional<T> string_to_cv_child(std::string_view s)
{
  return CV::from_string(s).and_then([](const auto& cv) -> std::optional<T> {
    return T(cv.code(), cv.accession());
  });
}

/******************************************************************************/
// clang doesn't support views::join yet :(
std::string join_(auto begin, auto end)
{
  std::string res(*begin);

  for (++begin; begin != end; ++begin) {
    res += "_" + *begin;
  }

  return res;
}

/******************************************************************************/
// Try to find a nested primitive node assuming that all fields will
// use the same type (i.e. homogeneous list).
std::shared_ptr<parquet::schema::PrimitiveNode>
find_homogeneous_primitive(std::shared_ptr<parquet::schema::Node> node)
{
  while (node != nullptr && node->is_group()) {
    auto group = std::static_pointer_cast<parquet::schema::GroupNode>(node);

    if (group->field_count() > 0) {
      node = group->field(0);
    } else {
      node = nullptr;
    }
  }

  if (node != nullptr && node->is_primitive()) {
    return std::static_pointer_cast<parquet::schema::PrimitiveNode>(node);
  } else {
    return nullptr;
  }
}

/******************************************************************************/
std::pair<Group::Field::Kind, std::optional<Util::Type>>
field_type_from_parquet(const std::shared_ptr<parquet::schema::GroupNode>& node)
{
  if (node->logical_type()->type() == parquet::LogicalType::Type::LIST) {
    // If the list is a homogeneous collection of scalars...
    if (auto prim = find_homogeneous_primitive(
            std::static_pointer_cast<parquet::schema::Node>(node));
        prim != nullptr) {
      return std::make_pair(Group::Field::Kind::List,
                            Util::type_from_parquet(*prim));
    } else {
      return std::make_pair(Group::Field::Kind::List, std::nullopt);
    }
  }

  return std::make_pair(Group::Field::Kind::Unknown, std::nullopt);
}

/******************************************************************************/
Group::Field::Field(std::string_view column_name,
                    index_type rel_index,
                    index_type abs_index)
    : rel_index_(rel_index)
    , abs_index_(abs_index)
    , schema_name_(column_name)
    , cv_type_()
    , cv_unit_()
    , type_()
{
}

/******************************************************************************/
Group::index_type Group::Field::relative_index() const { return rel_index_; }

/******************************************************************************/
Group::index_type Group::Field::absolute_index() const { return abs_index_; }

/******************************************************************************/
const std::string& Group::Field::name() const { return schema_name_; }

/******************************************************************************/
Group::Field::Kind Group::Field::kind() const { return kind_; }

/******************************************************************************/
const std::optional<Group::CVType>& Group::Field::cv_type() const
{
  return cv_type_;
}

/******************************************************************************/
const std::optional<Group::CVUnit>& Group::Field::cv_unit() const
{
  return cv_unit_;
}

/******************************************************************************/
const std::optional<Util::Type>& Group::Field::type() const { return type_; }

/******************************************************************************/
void Group::Field::type(Util::Type type) { type_ = type; }

/******************************************************************************/
std::optional<Util::Numpress::Type> Group::Field::possibly_numpress() const
{
  return Util::Numpress::type_from_column_name(schema_name_);
}

/******************************************************************************/
Group::Group(const parquet::schema::GroupNode& node, const Schema::File& file)
    : name_("root")
    , is_root_(true)
    , index_(0)
    , fields_()
{
  // This is the root group so only collect non-group top-level columns.
  make_fields(node, file, 0);
}

/******************************************************************************/
Group::Group(const parquet::schema::GroupNode& node,
             const Schema::File& file,
             index_type index,
             index_type offset)
    : name_(node.name())
    , is_root_(false)
    , index_(index)
    , fields_()
{
  make_fields(node, file, offset);
}

/******************************************************************************/
void Group::make_fields(const parquet::schema::GroupNode& node,
                        const Schema::File& file,
                        index_type offset)
{
  auto link = [&](std::shared_ptr<Field>& field) -> void {
    fields_[field->name()] = field;

    std::string col_path = path(*field);
    const auto it = std::ranges::find(file.columns(), col_path, &File::Column::path);

    if (it != file.columns().end()) {
      field->cv_type_ = it->accession.and_then(&string_to_cv_child<CVType>);
      field->cv_unit_ = it->unit.and_then(&string_to_cv_child<CVUnit>);
    }
  };

  for (index_type i : std::views::iota(0, node.field_count())) {
    auto child = node.field(i);

    if (child->is_primitive()) {
      std::shared_ptr<Field> field =
          std::make_shared<Field>(child->name(), i, offset + i);
      link(field);

      auto prim = std::static_pointer_cast<parquet::schema::PrimitiveNode>(child);
      field->kind_ = Field::Kind::Scalar;
      field->type_ = Util::type_from_parquet(*prim);
    } else if (child->is_group()) {
      auto grp = std::static_pointer_cast<parquet::schema::GroupNode>(child);

      if (grp->field_count() == 1) {
        std::shared_ptr<Field> field =
            std::make_shared<Field>(grp->name(), i, offset + i);
        link(field);

        auto grp_type = field_type_from_parquet(grp);
        field->kind_ = grp_type.first;
        field->type_ = grp_type.second;

        if (field->name() == "parameters") {
          field->kind_ = Group::Field::Kind::Params;
        }
      } else if (grp->field_count() > 1) {
        // A group with several children is a STRUCT.  Only the single-child
        // shape -- Parquet's LIST wrapper -- was registered before, so a struct
        // column was skipped silently and any array index entry naming it
        // failed to resolve at all.  The grid encoding's `<array>_grid` is such
        // a column; it is handed to the caller whole, as the struct it is,
        // because its three children are only meaningful together.
        std::shared_ptr<Field> field =
            std::make_shared<Field>(grp->name(), i, offset + i);
        link(field);
        field->kind_ = Group::Field::Kind::Struct;
      }
    }
  }
}

/******************************************************************************/
const std::string& Group::name() const { return name_; }

/******************************************************************************/
Group::index_type Group::index() const { return index_; }

/******************************************************************************/
std::optional<std::shared_ptr<const Group::Field>>
Group::field(std::string_view name) const
{
  auto it = fields_.find(std::string{name});

  if (it != fields_.end()) {
    return it->second;
  } else {
    return {};
  }
}

/******************************************************************************/
std::optional<std::shared_ptr<const Group::Field>>
Group::field(const CVType&& cvt) const
{
  const auto fields = fields_ | std::views::values;

  auto it =
      std::ranges::find_if(fields, [cvt](const std::shared_ptr<Field>& f) -> bool {
        return f->cv_type().has_value() && (f->cv_type().value() == cvt);
      });

  if (it == std::ranges::end(fields)) {
    return {};
  } else {
    return *it;
  }
}

/******************************************************************************/
const Group::field_map_t& Group::fields() const { return fields_; }

/******************************************************************************/
bool Group::is_root() const { return is_root_; }

/******************************************************************************/
std::string Group::path(const Field& field) const
{
  if (is_root_) {
    return field.name();
  } else {
    return name() + "." + field.name();
  }
}

} // namespace MzPeak::Schema
