/**
 * @file std_string_bytes_caster.h
 * @brief Extended std::string type caster that also accepts Python bytes
 *
 * This overrides nanobind's default std::string caster to additionally
 * accept bytes objects so that both bytes and str are accepted for
 * string parameters.
 */

#pragma once

#include <nanobind/nanobind.h>
#include <string>

// We must NOT include <nanobind/stl/string.h> when using this caster.
// Instead, we provide our own specialization.

namespace nanobind {
namespace detail {

template <>
struct type_caster<std::string> {
    NB_TYPE_CASTER(std::string, const_name("str"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
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
            value = std::string(data, (size_t)size);
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
            value = std::string(data, (size_t)size);
            return true;
        }

        return false;
    }

    static handle from_cpp(const std::string& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return PyUnicode_FromStringAndSize(src.c_str(), src.size());
    }

    static handle from_cpp(std::string&& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return PyUnicode_FromStringAndSize(src.c_str(), src.size());
    }
};

}  // namespace detail
}  // namespace nanobind
