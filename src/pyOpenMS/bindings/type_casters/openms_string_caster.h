/**
 * @file string.h
 * @brief nanobind type caster for OpenMS::String
 *
 * Provides bidirectional conversion between Python str/bytes and OpenMS::String.
 */

#pragma once

#include <nanobind/nanobind.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace nanobind {
namespace detail {

/**
 * Type caster for OpenMS::String <-> Python str
 *
 * This caster enables automatic conversion between Python strings (str)
 * and OpenMS::String in both directions.
 *
 * Python -> C++:
 *   - str: decoded as UTF-8 into OpenMS::String
 *   - bytes: used directly as OpenMS::String content
 *
 * C++ -> Python:
 *   - OpenMS::String: converted to Python str (assuming UTF-8)
 */
template <>
struct type_caster<OpenMS::String> {
public:
    // Python type hint
    NB_TYPE_CASTER(OpenMS::String, const_name("str"))

    /**
     * Convert Python object to OpenMS::String
     */
    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        // Handle None
        if (src.is_none()) {
            return false;
        }

        // Try str first
        if (PyUnicode_Check(src.ptr())) {
            Py_ssize_t size;
            const char* data = PyUnicode_AsUTF8AndSize(src.ptr(), &size);
            if (!data) {
                PyErr_Clear();
                return false;
            }
            value = OpenMS::String(data, size);
            return true;
        }

        // Try bytes
        if (PyBytes_Check(src.ptr())) {
            char* data;
            Py_ssize_t size;
            if (PyBytes_AsStringAndSize(src.ptr(), &data, &size) == -1) {
                PyErr_Clear();
                return false;
            }
            value = OpenMS::String(data, size);
            return true;
        }

        // Try objects with __str__ method
        PyObject* str_obj = PyObject_Str(src.ptr());
        if (str_obj) {
            Py_ssize_t size;
            const char* data = PyUnicode_AsUTF8AndSize(str_obj, &size);
            if (data) {
                value = OpenMS::String(data, size);
                Py_DECREF(str_obj);
                return true;
            }
            Py_DECREF(str_obj);
            PyErr_Clear();
        }
        PyErr_Clear();

        return false;
    }

    /**
     * Convert OpenMS::String to Python object
     */
    static handle from_cpp(const OpenMS::String& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return PyUnicode_FromStringAndSize(src.c_str(), src.size());
    }

    /**
     * Convert OpenMS::String reference to Python object
     */
    static handle from_cpp(OpenMS::String& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::String&>(src), policy, cleanup);
    }

    /**
     * Convert OpenMS::String rvalue to Python object
     */
    static handle from_cpp(OpenMS::String&& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::String&>(src), policy, cleanup);
    }
};

/**
 * Type caster for const OpenMS::String&
 * Delegates to the main caster
 */
template <>
struct type_caster<const OpenMS::String&> : type_caster<OpenMS::String> {};

/**
 * Type caster for OpenMS::String*
 */
template <>
struct type_caster<OpenMS::String*> {
    NB_TYPE_CASTER(OpenMS::String*, const_name("str | None"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        if (src.is_none()) {
            value = nullptr;
            return true;
        }

        type_caster<OpenMS::String> inner;
        if (!inner.from_python(src, flags, cleanup)) {
            return false;
        }
        // Note: This requires proper lifetime management
        // The caller must ensure the string outlives its use
        temp_value_ = std::move(inner.value);
        value = &temp_value_;
        return true;
    }

    static handle from_cpp(OpenMS::String* src, rv_policy policy, cleanup_list* cleanup) noexcept {
        if (!src) {
            return none().release();
        }
        return type_caster<OpenMS::String>::from_cpp(*src, policy, cleanup);
    }

private:
    OpenMS::String temp_value_;
};

}  // namespace detail
}  // namespace nanobind
