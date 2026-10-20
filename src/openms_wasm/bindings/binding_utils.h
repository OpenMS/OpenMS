// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// Helpers shared by the Embind binding files of the OpenMS JavaScript module.

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ParamValue.h>

#include "js_types.h"

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace OpenMS::WASM
{
  using emscripten::val;

  //------------------------------------------------------------------------------
  // Exceptions
  //
  // A C++ exception that leaves the module reaches JavaScript as an opaque
  // WebAssembly.Exception, without its message. Every bound function therefore goes
  // through Guard<>, which turns a C++ exception into a JavaScript Error carrying the
  // OpenMS exception name and message (e.g. "ParseError: ...").

  /// The message of the exception currently being handled (call from a catch block).
  inline std::string currentExceptionMessage()
  {
    try
    {
      throw;
    }
    catch (const Exception::BaseException& e)
    {
      return std::string(e.getName()) + ": " + e.what();
    }
    // Memory: a wasm32 module addresses at most 4 GB, a wasm64 one 16 GB (MAXIMUM_MEMORY, cmake/emscripten.cmake;
    // and its memory does not shrink again)
    catch (const std::bad_alloc&)
    {
#ifdef __wasm64__
      return "OutOfMemory: an allocation failed (the memory of this 64-bit WebAssembly module can grow to at most 16 GB)";
#else
      return "OutOfMemory: an allocation failed (a WebAssembly module can address at most 4 GB)";
#endif
    }
    catch (const std::length_error& e)
    {
#ifdef __wasm64__
      return std::string("OutOfMemory: ") + e.what() + ": a container would exceed what the module's memory can hold (16 GB)";
#else
      return std::string("OutOfMemory: ") + e.what() + ": a container would exceed what a 32-bit WebAssembly module can address";
#endif
    }
    catch (const std::exception& e)
    {
      return e.what();
    }
    catch (...)
    {
      return "unknown C++ exception";
    }
  }

  /// Throws a JavaScript Error with @p message. Never call it inside a C++ catch block:
  /// the C++ exception object would never be released.
  [[noreturn]] inline void throwJSError(const std::string& message)
  {
    val::global("Error").new_(message).throw_();
  }

  //------------------------------------------------------------------------------
  // Integers: numbers in JavaScript, and checked arguments
  //
  // Embind passes a 64-bit integer as BigInt. long and unsigned long (size_t, ptrdiff_t, ssize_t: Size,
  // SignedSize) are 64 bits in the wasm64 build and 32 bits in wasm32, so the guards below take and return
  // them as double in both builds: the JavaScript API and openms.d.ts are the same (number), and a result
  // is exact up to 2^53. Integer arguments (these and int, unsigned, short, char: Int, UInt, ...) are taken
  // as double as well, where Embind would truncate the number silently: an argument has to be an integer
  // within the range of its type (and of 2^53), else the call throws InvalidValue. (int64_t and uint64_t
  // are long long in both builds and BigInt in both.) The build checks that long and unsigned long cross nowhere else, e.g. as value of a val
  // (OPENMS_WASM_CHECK_INTEGERS, src/openms_wasm/CMakeLists.txt).

  /// Whether @p T is long or unsigned long, which Embind passes as BigInt in the wasm64 build
  template <typename T>
  inline constexpr bool isPointerSizedInteger = std::is_same_v<T, long> || std::is_same_v<T, unsigned long>;

  /// Whether @p T is an integer type whose arguments are checked: the ones Embind passes as number (in wasm64 after
  /// the conversion of long and unsigned long), i.e. not bool and not long long (BigInt in both builds). Embind
  /// truncates a number for them (1.5 to 1, -1 to 4294967295 for unsigned, 2^31 to -2^31 for int, NaN to 0).
  template <typename T>
  inline constexpr bool isCheckedInteger = isPointerSizedInteger<T> || std::is_same_v<T, int> || std::is_same_v<T, unsigned int>
    || std::is_same_v<T, short> || std::is_same_v<T, unsigned short> || std::is_same_v<T, char> || std::is_same_v<T, signed char>
    || std::is_same_v<T, unsigned char>;

  template <typename T>
  using Plain = std::remove_cv_t<std::remove_reference_t<T>>;

  /// Whether an argument of type @p T is taken as number and checked (integers but long long and bool)
  template <typename T>
  inline constexpr bool checksArgument = isCheckedInteger<Plain<T>>;

  /// Whether a result of type @p T is returned as double (long and unsigned long)
  template <typename T>
  inline constexpr bool convertsToNumber = isPointerSizedInteger<Plain<T>>;

  /// The type Embind sees for an argument of type @p T: double for the integers checked
  template <typename T>
  using JSArgument = std::conditional_t<checksArgument<T>, double, T>;

  /// The type Embind sees for a result of type @p T: double for long and unsigned long
  template <typename T>
  using JSResult = std::conditional_t<convertsToNumber<T>, double, T>;

  /// The integer of type @p I for the number @p value, an argument of a bound function: InvalidValue unless it
  /// is an integer within the range of @p I and of 2^53 (the safe integers). (The message names no position: a
  /// function bound as method takes the object as first C++ argument, which JavaScript does not count. Sizes and
  /// indices, size_t, are "a size or an index", other integers "an integer argument", unless @p what names it.)
  template <typename I>
  I integerFromJS(double value, const char* what = nullptr)
  {
    constexpr double safe = 9007199254740991.0; // 2^53 - 1
    constexpr double lowest = std::is_signed_v<I> ? std::max(-safe, static_cast<double>(std::numeric_limits<I>::min())) : 0.0;
    constexpr double highest = std::min(safe, static_cast<double>(std::numeric_limits<I>::max()));
    if (!(value >= lowest && value <= highest && std::trunc(value) == value)) // NaN fails
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        std::string(what ? what : std::is_same_v<I, unsigned long> ? "a size or an index" : "an integer argument") + " has to be an integer from "
        + std::to_string(static_cast<long long>(lowest)) + " to " + std::to_string(static_cast<unsigned long long>(highest)),
        val(value).call<std::string>("toString")); // NaN, Infinity, 1e+21 as JavaScript writes them
    }
    return static_cast<I>(value);
  }

  /// An argument of type @p A, from the value Embind passed (of type JSArgument<A>)
  template <typename A, typename J>
  decltype(auto) argumentFromJS(J&& value)
  {
    if constexpr (checksArgument<A>)
    {
      static_assert(!std::is_lvalue_reference_v<A> || std::is_const_v<std::remove_reference_t<A>>,
                    "a non-const reference to an integer cannot be bound");
      return integerFromJS<Plain<A>>(value);
    }
    else
    {
      return static_cast<J&&>(value);
    }
  }

  /// Calls @p fn with @p args (of the types JSArgument<A>...) as @p A..., and returns its result (of type @p R) as
  /// JSResult<R>
  template <typename R, typename... A, typename Fn, typename... J>
  JSResult<R> callFromJS(Fn&& fn, J&&... args)
  {
    if constexpr (std::is_void_v<R>)
    {
      fn(argumentFromJS<A>(static_cast<J&&>(args))...);
    }
    else if constexpr (convertsToNumber<R>)
    {
      const Plain<R> result = fn(argumentFromJS<A>(static_cast<J&&>(args))...);
#ifdef __wasm64__
      // A number is exact up to 2^53 (wasm32: always). Beyond, e.g. a sentinel such as std::string::npos (which is
      // 2^32 - 1 in wasm32), the binding has to map it to a number of its own (-1): it is an error, not a rounded value.
      if (result > static_cast<Plain<R>>(9007199254740991) || (std::is_signed_v<Plain<R>> && result < -static_cast<Plain<R>>(9007199254740991)))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "a size or an index returned to JavaScript has to be at most 2^53 - 1 (bindings/binding_utils.h)", std::to_string(result));
      }
#endif
      return static_cast<double>(result);
    }
    else
    {
      return fn(argumentFromJS<A>(static_cast<J&&>(args))...);
    }
  }

  /// Whether a call with arguments of the types @p A... checks one (and so can throw)
  template <typename... A>
  inline constexpr bool checksArguments = (checksArgument<A> || ...);

  /// Wraps the function (or member function) @p F so that a C++ exception thrown by it
  /// becomes a JavaScript Error. Guard<F>::call has the signature Embind expects: a member
  /// function becomes a free function taking the object as first argument. Integers of
  /// pointer size are numbers, and integer arguments are checked (see above).
  template <auto F, typename = decltype(F)>
  struct Guard;

#define OPENMS_WASM_GUARD_CALL(...)                     \
  std::string message;                                  \
  try                                                   \
  {                                                     \
    return __VA_ARGS__;                                 \
  }                                                     \
  catch (...)                                           \
  {                                                     \
    message = currentExceptionMessage();                \
  }                                                     \
  throwJSError(message);

/// The call of a guarded function: @p CALL (a lambda) with the arguments args... of the types JSArgument<A>...
#define OPENMS_WASM_CALL(CALL) callFromJS<R, A...>(CALL, std::forward<JSArgument<A>>(args)...)

/// The same for a noexcept function: unguarded, unless an argument is checked
#define OPENMS_WASM_NOEXCEPT_CALL(CALL)                 \
  if constexpr (checksArguments<A...>)                  \
  {                                                     \
    OPENMS_WASM_GUARD_CALL(OPENMS_WASM_CALL(CALL))      \
  }                                                     \
  else                                                  \
  {                                                     \
    return OPENMS_WASM_CALL(CALL);                      \
  }

/// A call of the free function F
#define OPENMS_WASM_FREE [](auto&&... a) -> decltype(auto) { return F(std::forward<decltype(a)>(a)...); }
/// A call of the member function F of self
#define OPENMS_WASM_MEMBER [&self](auto&&... a) -> decltype(auto) { return (self.*F)(std::forward<decltype(a)>(a)...); }

  template <auto F, typename R, typename... A>
  struct Guard<F, R (*)(A...)>
  {
    static JSResult<R> call(JSArgument<A>... args) { OPENMS_WASM_GUARD_CALL(OPENMS_WASM_CALL(OPENMS_WASM_FREE)) }
  };

  template <auto F, typename R, typename... A>
  struct Guard<F, R (*)(A...) noexcept>
  {
    static JSResult<R> call(JSArgument<A>... args) { OPENMS_WASM_NOEXCEPT_CALL(OPENMS_WASM_FREE) }
  };

  template <auto F, typename R, typename C, typename... A>
  struct Guard<F, R (C::*)(A...)>
  {
    static JSResult<R> call(C& self, JSArgument<A>... args) { OPENMS_WASM_GUARD_CALL(OPENMS_WASM_CALL(OPENMS_WASM_MEMBER)) }
  };

  template <auto F, typename R, typename C, typename... A>
  struct Guard<F, R (C::*)(A...) const>
  {
    static JSResult<R> call(const C& self, JSArgument<A>... args) { OPENMS_WASM_GUARD_CALL(OPENMS_WASM_CALL(OPENMS_WASM_MEMBER)) }
  };

  template <auto F, typename R, typename C, typename... A>
  struct Guard<F, R (C::*)(A...) noexcept>
  {
    static JSResult<R> call(C& self, JSArgument<A>... args) { OPENMS_WASM_NOEXCEPT_CALL(OPENMS_WASM_MEMBER) }
  };

  template <auto F, typename R, typename C, typename... A>
  struct Guard<F, R (C::*)(A...) const noexcept>
  {
    static JSResult<R> call(const C& self, JSArgument<A>... args) { OPENMS_WASM_NOEXCEPT_CALL(OPENMS_WASM_MEMBER) }
  };

  /// Guard<> for a member function @p F of class @p Self, which may be declared in a (public)
  /// base class: Embind needs the bound class itself as type of the object argument.
  template <typename Self, auto F, typename = decltype(F)>
  struct MethodGuard;

  template <typename Self, auto F, typename R, typename C, typename... A>
  struct MethodGuard<Self, F, R (C::*)(A...)>
  {
    static JSResult<R> call(Self& self, JSArgument<A>... args) { OPENMS_WASM_GUARD_CALL(OPENMS_WASM_CALL(OPENMS_WASM_MEMBER)) }
  };

  template <typename Self, auto F, typename R, typename C, typename... A>
  struct MethodGuard<Self, F, R (C::*)(A...) const>
  {
    static JSResult<R> call(const Self& self, JSArgument<A>... args) { OPENMS_WASM_GUARD_CALL(OPENMS_WASM_CALL(OPENMS_WASM_MEMBER)) }
  };

  template <typename Self, auto F, typename R, typename C, typename... A>
  struct MethodGuard<Self, F, R (C::*)(A...) noexcept>
  {
    static JSResult<R> call(Self& self, JSArgument<A>... args) { OPENMS_WASM_NOEXCEPT_CALL(OPENMS_WASM_MEMBER) }
  };

  template <typename Self, auto F, typename R, typename C, typename... A>
  struct MethodGuard<Self, F, R (C::*)(A...) const noexcept>
  {
    static JSResult<R> call(const Self& self, JSArgument<A>... args) { OPENMS_WASM_NOEXCEPT_CALL(OPENMS_WASM_MEMBER) }
  };

  /// A constructor of @p T taking @p A... (for class_::constructor()), guarded like the functions above: a
  /// C++ exception becomes a JavaScript Error, integers of pointer size are numbers, and integer arguments are checked.
  template <typename T, typename... A>
  struct ConstructorGuard
  {
    using R = T*;
    static T* call(JSArgument<A>... args)
    {
      OPENMS_WASM_GUARD_CALL(OPENMS_WASM_CALL([](auto&&... a) { return new T(std::forward<decltype(a)>(a)...); }))
    }
  };

#undef OPENMS_WASM_MEMBER
#undef OPENMS_WASM_FREE
#undef OPENMS_WASM_NOEXCEPT_CALL
#undef OPENMS_WASM_CALL
#undef OPENMS_WASM_GUARD_CALL

/// The guarded version of a free function (taking the object as first argument when bound
/// as method), for class_::function(), class_::class_function() and function().
#define OMS_FN(f) (&::OpenMS::WASM::Guard<f>::call)

/// The guarded version of member function @p f, bound as method of class @p Self.
#define OMS_METHOD(Self, f) (&::OpenMS::WASM::MethodGuard<Self, f>::call)

/// The guarded constructor of class @p T taking arguments of the types ..., for class_::constructor():
/// class_<T>(...).constructor(OMS_CONSTRUCTOR(T, Size)). (Unlike constructor<Size>(), it takes the size as number.)
#define OMS_CONSTRUCTOR(T, ...) (&::OpenMS::WASM::ConstructorGuard<T, __VA_ARGS__>::call)

  //------------------------------------------------------------------------------
  // Arrays

  template <typename T> struct TypedArrayName;
  template <> struct TypedArrayName<double> { static constexpr const char* value = "Float64Array"; };
  template <> struct TypedArrayName<float> { static constexpr const char* value = "Float32Array"; };
  template <> struct TypedArrayName<int32_t> { static constexpr const char* value = "Int32Array"; };
  template <> struct TypedArrayName<uint16_t> { static constexpr const char* value = "Uint16Array"; };
  template <> struct TypedArrayName<int16_t> { static constexpr const char* value = "Int16Array"; };
  template <> struct TypedArrayName<uint32_t> { static constexpr const char* value = "Uint32Array"; };
  template <> struct TypedArrayName<int8_t> { static constexpr const char* value = "Int8Array"; };
  template <> struct TypedArrayName<uint8_t> { static constexpr const char* value = "Uint8Array"; };

  /// A JavaScript typed array holding a copy of @p n values at @p data.
  template <typename T>
  val toTypedArray(const T* data, size_t n)
  {
    // the length as a number: a size_t of wasm64 reaches JavaScript as BigInt, which the constructor rejects
    val result = val::global(TypedArrayName<T>::value).new_(static_cast<double>(n));
    if (n > 0)
    {
      result.call<void>("set", val(emscripten::typed_memory_view(n, data)));
    }
    return result;
  }

  template <typename T>
  val toTypedArray(const std::vector<T>& v)
  {
    return toTypedArray(v.data(), v.size());
  }

  /// A JavaScript typed array of element type @p T with @p values (converted from their type).
  template <typename T, typename Range, typename Projection>
  val toTypedArray(const Range& range, Projection project)
  {
    std::vector<T> values;
    values.reserve(range.size());
    for (const auto& element : range)
    {
      values.push_back(static_cast<T>(project(element)));
    }
    return toTypedArray(values);
  }

  /// A JavaScript array of @p v's elements (converted by Embind: numbers, strings, value
  /// objects, or copies of bound classes).
  template <typename Range>
  val toArray(const Range& v)
  {
    val result = val::array();
    for (const auto& element : v)
    {
      result.call<void>("push", element);
    }
    return result;
  }

  /// A JavaScript array of @p v's elements, each converted by @p project.
  template <typename Range, typename Projection>
  val toArray(const Range& v, Projection project)
  {
    val result = val::array();
    for (const auto& element : v)
    {
      result.call<void>("push", project(element));
    }
    return result;
  }

  /// The length of a JavaScript array (or typed array, string): 0 for a value without one. (Read as number: a
  /// size_t would be a BigInt in wasm64, and an unsigned an unchecked conversion: OPENMS_WASM_CHECK_INTEGERS.)
  inline std::size_t lengthOf(const val& v)
  {
    const double length = v["length"].as<double>();
    return length > 0 ? static_cast<std::size_t>(length) : 0; // NaN: 0
  }

  /// Numbers of a JavaScript array or typed array (none for a value without length), as
  /// emscripten::convertJSArrayToNumberVector() (which reads the length as size_t). Integers (@p T int, unsigned,
  /// ...) have to be integers within the range of @p T, else InvalidValue (as arguments: integerFromJS()).
  template <typename T>
  std::vector<T> numbersFromJS(const val& v)
  {
    const std::size_t n = lengthOf(v);
    if constexpr (std::is_integral_v<T>)
    {
      const std::vector<double> numbers = numbersFromJS<double>(v);
      std::vector<T> result;
      result.reserve(n);
      for (const double number : numbers) result.push_back(integerFromJS<T>(number, "an element of a list of integers"));
      return result;
    }
    else
    {
      std::vector<T> result(n);
      if (n > 0)
      {
        val(emscripten::typed_memory_view(result.size(), result.data())).call<void>("set", v);
      }
      return result;
    }
  }

  /// The elements of a JavaScript array as @p T, e.g. strings or objects of a bound class (as
  /// emscripten::vecFromJSArray(), which reads the length as unsigned: an unchecked conversion).
  template <typename T>
  std::vector<T> vectorFromJS(const val& v)
  {
    const std::size_t n = lengthOf(v);
    std::vector<T> result;
    result.reserve(n);
    for (std::size_t i = 0; i < n; ++i) result.push_back(v[static_cast<double>(i)].as<T>());
    return result;
  }

  /// Strings of a JavaScript array.
  inline std::vector<std::string> stringsFromJS(const val& v)
  {
    return vectorFromJS<std::string>(v);
  }

  /// Whether @p v is undefined or null (an omitted optional argument).
  inline bool isNullish(const val& v)
  {
    return v.isUndefined() || v.isNull();
  }

  //------------------------------------------------------------------------------
  // Values of Param and DataValue (meta values)

  /// A Param or DataValue as JavaScript value: number, string, array, or null (empty).
  /// @p as_bool: a "true"/"false" string becomes a boolean (boolean parameters).
  val toJS(const ParamValue& v, bool as_bool = false);
  val toJS(const DataValue& v);

  /// A JavaScript value as ParamValue. @p type is the type of the value it replaces
  /// (EMPTY_VALUE: inferred from the JavaScript value; integral numbers become ints).
  ParamValue paramValueFromJS(const val& v, ParamValue::ValueType type = ParamValue::EMPTY_VALUE);
  DataValue dataValueFromJS(const val& v);

  /// All entries of @p p as {name: value} object.
  JS::ParamValues paramToObject(const Param& p);

  /// Sets the values in @p values (an object {name: value}); names must exist in @p p.
  void updateParamFromObject(Param& p, const JS::ParamValues& values);

  //------------------------------------------------------------------------------
  // Meta values (MetaInfoInterface), bound as methods of the classes that have them

  /// Meta values of @p object as {name: value} object.
  template <typename T>
  JS::MetaValues metaValuesToObject(const T& object)
  {
    val result = val::object();
    std::vector<std::string> keys;
    object.getKeys(keys);
    for (const auto& key : keys)
    {
      result.set(key, toJS(object.getMetaValue(key)));
    }
    return JS::MetaValues(result);
  }

  /// The meta value @p name, or undefined.
  template <typename T>
  JS::MetaValue metaValue(const T& object, const std::string& name)
  {
    return JS::MetaValue(object.metaValueExists(name) ? toJS(object.getMetaValue(name)) : val::undefined());
  }

  template <typename T>
  void setMetaValue(T& object, const std::string& name, const JS::MetaValueInput& value)
  {
    object.setMetaValue(name, dataValueFromJS(value));
  }

  template <typename T>
  bool metaValueExists(const T& object, const std::string& name)
  {
    return object.metaValueExists(name);
  }

  template <typename T>
  JS::StringArray metaValueKeys(const T& object)
  {
    std::vector<std::string> keys;
    object.getKeys(keys);
    return JS::StringArray(toArray(keys));
  }
} // namespace OpenMS::WASM
