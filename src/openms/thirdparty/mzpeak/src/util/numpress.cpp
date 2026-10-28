/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

// OpenMS: the MSNumpress copy that libOpenMS compiles (FORMAT/MSNUMPRESS), not the meson subproject
#include <OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h>

#include "mzpeak/exception.h"
#include "mzpeak/util/compat.h" // IWYU pragma: keep
#include "mzpeak/util/decoders.h"
#include "mzpeak/util/numpress.h"

namespace MzPeak::Util::Numpress {

/******************************************************************************/
using decoder =
    std::move_only_function<void(const std::vector<uint8_t>&, std::vector<double>&)>;

/******************************************************************************/
std::optional<Type> type_from_column_name(const std::string& name)
{
  if (name.contains("numpress_linear")) {
    return Util::Numpress::Linear;
  } else if (name.contains("numpress_slof")) {
    return Util::Numpress::SLOF;
  } else if (name.contains("numpress_pic")) {
    return Util::Numpress::PIC;
  } else {
    return {};
  }
}

/******************************************************************************/
std::size_t decoding_space_needed(std::size_t n, Type t)
{
  switch (t) {
  case Numpress::Linear:
    // Need C++26 for saturating_sub :(
    if (n <= 8) return 0;
    return (n - 8) * 2;
  case Numpress::SLOF:
    // FIXME: Could this be a typo in the numpress lib?
    // Need C++26 for saturating_sub :(
    if (n <= 8) return 0;
    return (n - 8) / 2;
  case Numpress::PIC:
    return n * 2;
  }

  std::unreachable();
}

/******************************************************************************/
std::shared_ptr<std::vector<double>>
from_arrow(const std::shared_ptr<arrow::Array>& src, decoder f)
{
  if (src->type_id() != arrow::Type::UINT8) {
    std::string msg("numpress decoding requested but source array is not uint8");
    throw InvalidFormatError(msg);
  }

  std::vector<uint8_t> bytes;
  bytes.reserve(src->length());

  Decoders::Scalar<uint8_t> decoder;
  decoder.decode(src, bytes);

  auto values = std::make_shared<std::vector<double>>();
  f(bytes, *values);

  return values;
}

/******************************************************************************/
void decode_linear(const std::vector<uint8_t>& input, std::vector<double>& output)
{
  try {
    ms::numpress::MSNumpress::decodeLinear(input, output);
  } catch (const char* msg) {
    throw InvalidFormatError(msg);
  }
}

/******************************************************************************/
std::shared_ptr<std::vector<double>>
decode_linear(const std::shared_ptr<arrow::Array>& src)
{
  return from_arrow(src, [](auto& i, auto& o) -> void { decode_linear(i, o); });
}

/******************************************************************************/
void decode_slof(const std::vector<uint8_t>& input, std::vector<double>& output)
{
  try {
    ms::numpress::MSNumpress::decodeSlof(input, output);
  } catch (const char* msg) {
    throw InvalidFormatError(msg);
  }
}

/******************************************************************************/
std::shared_ptr<std::vector<double>>
decode_slof(const std::shared_ptr<arrow::Array>& src)
{
  return from_arrow(src, [](auto& i, auto& o) -> void { decode_slof(i, o); });
}

/******************************************************************************/
void decode_pic(const std::vector<uint8_t>& input, std::vector<double>& output)
{
  try {
    ms::numpress::MSNumpress::decodePic(input, output);
  } catch (const char* msg) {
    throw InvalidFormatError(msg);
  }
}

/******************************************************************************/
std::shared_ptr<std::vector<double>>
decode_pic(const std::shared_ptr<arrow::Array>& src)
{
  return from_arrow(src, [](auto& i, auto& o) -> void { decode_pic(i, o); });
}

} // namespace MzPeak::Util::Numpress
