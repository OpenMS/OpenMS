/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/decoders.h"
#include "mzpeak/util/numpress.h"

namespace MzPeak::Util::Decoders {

/******************************************************************************/
std::size_t guess_array_length(const Schema::Column& column,
                               const std::shared_ptr<arrow::Array>& ary)
{
  std::optional<Numpress::Type> numpress = column.second->possibly_numpress();

  auto count =
      [&numpress](const std::shared_ptr<arrow::Array>& nums) -> std::size_t {
    if (numpress.has_value()) {
      return Numpress::decoding_space_needed(nums->length(), numpress.value());
    } else {
      return static_cast<std::size_t>(nums->length());
    }
  };

  auto for_list = [&]<typename T>(const std::shared_ptr<T>& list) -> std::size_t {
    std::size_t size = {};

    for (int64_t index : std::views::iota(0, list->length())) {
      if (list->IsValid(index)) {
        size += count(list->value_slice(index));
      }
    }

    return size;
  };

  auto type = ary->type_id();

  if (type == arrow::Type::LIST) {
    return for_list(std::static_pointer_cast<arrow::ListArray>(ary));
  } else if (type == arrow::Type::FIXED_SIZE_LIST) {
    return for_list(std::static_pointer_cast<arrow::FixedSizeListArray>(ary));
  } else if (type == arrow::Type::LARGE_LIST) {
    return for_list(std::static_pointer_cast<arrow::LargeListArray>(ary));
  } else if (type == arrow::Type::LIST_VIEW) {
    return for_list(std::static_pointer_cast<arrow::ListViewArray>(ary));
  } else if (type == arrow::Type::LARGE_LIST_VIEW) {
    return for_list(std::static_pointer_cast<arrow::LargeListViewArray>(ary));
  } else {
    return count(ary);
  }
}

} // namespace MzPeak::Util::Decoders
