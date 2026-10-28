/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

// Compatibility with compilers that don't fully support C++23.

/// std::move_only_function
#if defined(_LIBCPP_VERSION) && !defined(__cpp_lib_move_only_function)
#include <boost/compat/move_only_function.hpp>
namespace std {
template <class S> using move_only_function = boost::compat::move_only_function<S>;
} // namespace std
#endif

/******************************************************************************/
/// libzip takes file and member names as UTF-8 `const char*` on every platform,
/// while `path::c_str()` is `wchar_t*` on Windows and `path::string()` is the
/// ANSI code page there. One conversion for both uses; forward slashes are
/// what archive members need and what Win32 accepts for files.
#include <filesystem>
#include <string>

namespace MzPeak::Util {
inline std::string narrow(const std::filesystem::path& path)
{
  const std::u8string u8 = path.generic_u8string();
  return std::string(u8.begin(), u8.end());
}
} // namespace MzPeak::Util

/******************************************************************************/
// Helper to produce useful messages with static_assert.
template <typename...> inline constexpr bool false_type = false;

// #endif
