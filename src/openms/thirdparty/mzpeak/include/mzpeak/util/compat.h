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
/// A non-owning reference to a callable, as boost::compat::function_ref (Boost >= 1.86) and C++26's
/// std::function_ref: query.h's evaluation callbacks, so that the library builds with older Boost releases
/// (e.g. 1.83 of Ubuntu 24.04). The callable must outlive the reference (as a temporary argument does).
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

namespace MzPeak::Util {
template <class Signature> class function_ref;

template <class R, class... Args> class function_ref<R(Args...)>
{
public:
  template <class F>
    requires(!std::is_same_v<std::remove_cvref_t<F>, function_ref> && std::is_object_v<std::remove_reference_t<F>>
             && std::is_invocable_r_v<R, std::remove_reference_t<F>&, Args...>)
  function_ref(F&& f) noexcept
      : object_(const_cast<void*>(static_cast<const void*>(std::addressof(f)))),
        call_([](void* object, Args... args) -> R {
          return std::invoke_r<R>(*static_cast<std::remove_reference_t<F>*>(object), std::forward<Args>(args)...);
        })
  {
  }

  R operator()(Args... args) const { return call_(object_, std::forward<Args>(args)...); }

private:
  void* object_;
  R (*call_)(void*, Args...);
};
} // namespace MzPeak::Util

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
