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

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <cstdint>
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

  /// Wraps the function (or member function) @p F so that a C++ exception thrown by it
  /// becomes a JavaScript Error. Guard<F>::call has the signature Embind expects: a member
  /// function becomes a free function taking the object as first argument.
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

  template <auto F, typename R, typename... A>
  struct Guard<F, R (*)(A...)>
  {
    static R call(A... args) { OPENMS_WASM_GUARD_CALL(F(std::forward<A>(args)...)) }
  };

  template <auto F, typename R, typename... A>
  struct Guard<F, R (*)(A...) noexcept>
  {
    static R call(A... args) { return F(std::forward<A>(args)...); }
  };

  template <auto F, typename R, typename C, typename... A>
  struct Guard<F, R (C::*)(A...)>
  {
    static R call(C& self, A... args) { OPENMS_WASM_GUARD_CALL((self.*F)(std::forward<A>(args)...)) }
  };

  template <auto F, typename R, typename C, typename... A>
  struct Guard<F, R (C::*)(A...) const>
  {
    static R call(const C& self, A... args) { OPENMS_WASM_GUARD_CALL((self.*F)(std::forward<A>(args)...)) }
  };

  template <auto F, typename R, typename C, typename... A>
  struct Guard<F, R (C::*)(A...) noexcept>
  {
    static R call(C& self, A... args) { return (self.*F)(std::forward<A>(args)...); }
  };

  template <auto F, typename R, typename C, typename... A>
  struct Guard<F, R (C::*)(A...) const noexcept>
  {
    static R call(const C& self, A... args) { return (self.*F)(std::forward<A>(args)...); }
  };

  /// Guard<> for a member function @p F of class @p Self, which may be declared in a (public)
  /// base class: Embind needs the bound class itself as type of the object argument.
  template <typename Self, auto F, typename = decltype(F)>
  struct MethodGuard;

  template <typename Self, auto F, typename R, typename C, typename... A>
  struct MethodGuard<Self, F, R (C::*)(A...)>
  {
    static R call(Self& self, A... args) { OPENMS_WASM_GUARD_CALL((self.*F)(std::forward<A>(args)...)) }
  };

  template <typename Self, auto F, typename R, typename C, typename... A>
  struct MethodGuard<Self, F, R (C::*)(A...) const>
  {
    static R call(const Self& self, A... args) { OPENMS_WASM_GUARD_CALL((self.*F)(std::forward<A>(args)...)) }
  };

  template <typename Self, auto F, typename R, typename C, typename... A>
  struct MethodGuard<Self, F, R (C::*)(A...) noexcept>
  {
    static R call(Self& self, A... args) { return (self.*F)(std::forward<A>(args)...); }
  };

  template <typename Self, auto F, typename R, typename C, typename... A>
  struct MethodGuard<Self, F, R (C::*)(A...) const noexcept>
  {
    static R call(const Self& self, A... args) { return (self.*F)(std::forward<A>(args)...); }
  };

#undef OPENMS_WASM_GUARD_CALL

/// The guarded version of a free function (taking the object as first argument when bound
/// as method), for class_::function(), class_::class_function() and function().
#define OMS_FN(f) (&::OpenMS::WASM::Guard<f>::call)

/// The guarded version of member function @p f, bound as method of class @p Self.
#define OMS_METHOD(Self, f) (&::OpenMS::WASM::MethodGuard<Self, f>::call)

  //------------------------------------------------------------------------------
  // Arrays

  template <typename T> struct TypedArrayName;
  template <> struct TypedArrayName<double> { static constexpr const char* value = "Float64Array"; };
  template <> struct TypedArrayName<float> { static constexpr const char* value = "Float32Array"; };
  template <> struct TypedArrayName<int32_t> { static constexpr const char* value = "Int32Array"; };
  template <> struct TypedArrayName<uint32_t> { static constexpr const char* value = "Uint32Array"; };
  template <> struct TypedArrayName<int8_t> { static constexpr const char* value = "Int8Array"; };
  template <> struct TypedArrayName<uint8_t> { static constexpr const char* value = "Uint8Array"; };

  /// A JavaScript typed array holding a copy of @p n values at @p data.
  template <typename T>
  val toTypedArray(const T* data, size_t n)
  {
    val result = val::global(TypedArrayName<T>::value).new_(n);
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

  /// Numbers of a JavaScript array or typed array.
  template <typename T>
  std::vector<T> numbersFromJS(const val& v)
  {
    return emscripten::convertJSArrayToNumberVector<T>(v);
  }

  /// Strings of a JavaScript array.
  inline std::vector<std::string> stringsFromJS(const val& v)
  {
    return emscripten::vecFromJSArray<std::string>(v);
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
  val paramToObject(const Param& p);

  /// Sets the values in @p values (an object {name: value}); names must exist in @p p.
  void updateParamFromObject(Param& p, const val& values);

  /// Meta values of a MetaInfoInterface as {name: value} object.
  template <typename T>
  val metaValuesToObject(const T& object)
  {
    val result = val::object();
    std::vector<std::string> keys;
    object.getKeys(keys);
    for (const auto& key : keys)
    {
      result.set(key, toJS(object.getMetaValue(key)));
    }
    return result;
  }
} // namespace OpenMS::WASM
