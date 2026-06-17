/**
 * @file string.h
 * @brief nanobind type caster for std::string
 *
 * Provides bidirectional conversion between Python str/bytes and std::string.
 */

#pragma once

#include <nanobind/nanobind.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

namespace nanobind {
namespace detail {

/**
 * Type caster for std::string <-> Python str
 *
 * This caster enables automatic conversion between Python strings (str)
 * and std::string in both directions.
 *
 * Python -> C++:
 *   - str: decoded as UTF-8 into std::string
 *   - bytes: used directly as std::string content
 *
 * C++ -> Python:
 *   - std::string: converted to Python str (assuming UTF-8)
 */
template <>
struct type_caster<std::string> {
public:
    // Python type hint
    NB_TYPE_CASTER(std::string, const_name("str"))

    /**
     * Convert Python object to std::string
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
            value = std::string(data, size);
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
            value = std::string(data, size);
            return true;
        }

        // Try objects with __str__ method
        PyObject* str_obj = PyObject_Str(src.ptr());
        if (str_obj) {
            Py_ssize_t size;
            const char* data = PyUnicode_AsUTF8AndSize(str_obj, &size);
            if (data) {
                value = std::string(data, size);
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
     * Convert std::string to Python object
     */
    static handle from_cpp(const std::string& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return PyUnicode_FromStringAndSize(src.c_str(), src.size());
    }

    /**
     * Convert std::string reference to Python object
     */
    static handle from_cpp(std::string& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const std::string&>(src), policy, cleanup);
    }

    /**
     * Convert std::string rvalue to Python object
     */
    static handle from_cpp(std::string&& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const std::string&>(src), policy, cleanup);
    }
};

/**
 * Type caster for const std::string&
 * Delegates to the main caster
 */
template <>
struct type_caster<const std::string&> : type_caster<std::string> {};

/**
 * Type caster for std::string*
 */
template <>
struct type_caster<std::string*> {
    NB_TYPE_CASTER(std::string*, const_name("str | None"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        if (src.is_none()) {
            value = nullptr;
            return true;
        }

        type_caster<std::string> inner;
        if (!inner.from_python(src, flags, cleanup)) {
            return false;
        }
        // Note: This requires proper lifetime management
        // The caller must ensure the string outlives its use
        temp_value_ = std::move(inner.value);
        value = &temp_value_;
        return true;
    }

    static handle from_cpp(std::string* src, rv_policy policy, cleanup_list* cleanup) noexcept {
        if (!src) {
            return none().release();
        }
        return type_caster<std::string>::from_cpp(*src, policy, cleanup);
    }

private:
    std::string temp_value_;
};

}  // namespace detail
}  // namespace nanobind
