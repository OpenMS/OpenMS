/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <utility>

#include "mzpeak/schema/psi/chunk_encoding.h"

namespace MzPeak::Schema::PSI {

/******************************************************************************/
ChunkEncoding::ChunkEncoding(const CV& cv)
    : val_(cv)
{
  if (cv.code() == "MS") {
    const auto& accession = cv.accession();

    if (accession == "1000576") {
      val_ = NoCompression;
    } else if (accession == "1003089") {
      val_ = Delta;
    } else if (accession == "1002312") {
      val_ = NumpressLinear;
    } else if (accession == "1003826") {
      val_ = Grid;
    }
  }
}

/******************************************************************************/
CV ChunkEncoding::to_cv() const
{
  return std::visit(
      [](auto&& v) -> CV {
        using T = std::decay_t<decltype(v)>;
        if constexpr (std::is_same_v<T, CV>) {
          return v;
        } else if constexpr (std::is_same_v<T, Type>) {
          switch (v) {
          case NoCompression:
            return CV("MS", "1000576");
          case Delta:
            return CV("MS", "1003089");
          case NumpressLinear:
            return CV("MS", "1002312");
          case Grid:
            return CV("MS", "1003826");
          }

          std::unreachable();
        }
      },
      val_);
}

/******************************************************************************/
std::optional<ChunkEncoding::Type> ChunkEncoding::type() const
{
  if (std::holds_alternative<Type>(val_)) {
    return std::get<Type>(val_);
  } else {
    return {};
  }
}

} // namespace MzPeak::Schema::PSI
