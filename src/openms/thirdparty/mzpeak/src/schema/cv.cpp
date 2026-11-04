/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include "mzpeak/schema/cv.h"

namespace MzPeak::Schema {

/******************************************************************************/
std::optional<CV> CV::from_string(std::string_view s)
{
  std::string_view::size_type sep_pos = s.find(':');

  if (sep_pos == std::string_view::npos) {
    return {};
  }

  return CV(s.substr(0, sep_pos), s.substr(sep_pos + 1));
}

/******************************************************************************/
std::string CV::to_string() const { return code_ + ":" + accession_; }

} // namespace MzPeak::Schema
