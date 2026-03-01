/**
 * @file stl_containers.h
 * @brief nanobind type casters for OpenMS STL container specializations
 *
 * Provides type casters for:
 *   - std::vector<OpenMS::String>
 *   - std::set<OpenMS::String>
 *   - std::map<OpenMS::String, ...>
 *   - OpenMS::StringList, IntList, DoubleList
 *   - CVTermMap (nested map structure)
 */

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/set.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/CVTerm.h>

#include "string.h"

namespace nanobind {
namespace detail {

/**
 * Type caster for std::vector<OpenMS::String>
 *
 * Optimized conversion between Python list[str] and C++ vector<String>
 */
template <>
struct type_caster<std::vector<OpenMS::String>> {
public:
    NB_TYPE_CASTER(std::vector<OpenMS::String>, const_name("list[str]"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        if (src.is_none()) {
            value.clear();
            return true;
        }

        if (!PySequence_Check(src.ptr())) {
            return false;
        }

        Py_ssize_t len = PySequence_Length(src.ptr());
        if (len < 0) {
            PyErr_Clear();
            return false;
        }

        value.clear();
        value.reserve(len);

        for (Py_ssize_t i = 0; i < len; ++i) {
            PyObject* item = PySequence_GetItem(src.ptr(), i);
            if (!item) {
                PyErr_Clear();
                return false;
            }

            if (PyUnicode_Check(item)) {
                Py_ssize_t size;
                const char* data = PyUnicode_AsUTF8AndSize(item, &size);
                if (data) {
                    value.push_back(OpenMS::String(data, size));
                } else {
                    Py_DECREF(item);
                    PyErr_Clear();
                    return false;
                }
            } else if (PyBytes_Check(item)) {
                char* data;
                Py_ssize_t size;
                if (PyBytes_AsStringAndSize(item, &data, &size) == 0) {
                    value.push_back(OpenMS::String(data, size));
                } else {
                    Py_DECREF(item);
                    PyErr_Clear();
                    return false;
                }
            } else {
                // Try __str__
                PyObject* str = PyObject_Str(item);
                if (str) {
                    Py_ssize_t size;
                    const char* data = PyUnicode_AsUTF8AndSize(str, &size);
                    if (data) {
                        value.push_back(OpenMS::String(data, size));
                    }
                    Py_DECREF(str);
                } else {
                    Py_DECREF(item);
                    PyErr_Clear();
                    return false;
                }
            }

            Py_DECREF(item);
        }

        return true;
    }

    static handle from_cpp(const std::vector<OpenMS::String>& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        PyObject* list = PyList_New(src.size());
        if (!list) return handle();

        for (size_t i = 0; i < src.size(); ++i) {
            PyObject* item = PyUnicode_FromStringAndSize(src[i].c_str(), src[i].size());
            if (!item) {
                Py_DECREF(list);
                return handle();
            }
            PyList_SET_ITEM(list, i, item);
        }

        return list;
    }

    static handle from_cpp(std::vector<OpenMS::String>& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const std::vector<OpenMS::String>&>(src), policy, cleanup);
    }

    static handle from_cpp(std::vector<OpenMS::String>&& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const std::vector<OpenMS::String>&>(src), policy, cleanup);
    }
};

/**
 * Type caster for std::set<OpenMS::String>
 */
template <>
struct type_caster<std::set<OpenMS::String>> {
public:
    NB_TYPE_CASTER(std::set<OpenMS::String>, const_name("set[str]"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        if (src.is_none()) {
            value.clear();
            return true;
        }

        // Accept set, frozenset, or any iterable
        PyObject* iter = PyObject_GetIter(src.ptr());
        if (!iter) {
            PyErr_Clear();
            return false;
        }

        value.clear();

        PyObject* item;
        while ((item = PyIter_Next(iter))) {
            if (PyUnicode_Check(item)) {
                Py_ssize_t size;
                const char* data = PyUnicode_AsUTF8AndSize(item, &size);
                if (data) {
                    value.insert(OpenMS::String(data, size));
                }
            } else if (PyBytes_Check(item)) {
                char* data;
                Py_ssize_t size;
                if (PyBytes_AsStringAndSize(item, &data, &size) == 0) {
                    value.insert(OpenMS::String(data, size));
                }
            } else {
                Py_DECREF(item);
                Py_DECREF(iter);
                return false;  // reject non-string items
            }
            Py_DECREF(item);
        }

        Py_DECREF(iter);

        if (PyErr_Occurred()) {
            PyErr_Clear();
            return false;
        }

        return true;
    }

    static handle from_cpp(const std::set<OpenMS::String>& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        PyObject* set = PySet_New(nullptr);
        if (!set) return handle();

        for (const auto& s : src) {
            PyObject* item = PyUnicode_FromStringAndSize(s.c_str(), s.size());
            if (!item || PySet_Add(set, item) == -1) {
                Py_XDECREF(item);
                Py_DECREF(set);
                return handle();
            }
            Py_DECREF(item);
        }

        return set;
    }

    static handle from_cpp(std::set<OpenMS::String>& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const std::set<OpenMS::String>&>(src), policy, cleanup);
    }

    static handle from_cpp(std::set<OpenMS::String>&& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const std::set<OpenMS::String>&>(src), policy, cleanup);
    }
};

/**
 * Type caster for std::map<OpenMS::String, V>
 *
 * This is a template that handles maps with OpenMS::String keys.
 */
template <typename V>
struct type_caster<std::map<OpenMS::String, V>> {
public:
    using MapType = std::map<OpenMS::String, V>;
    NB_TYPE_CASTER(MapType, const_name("dict[str, ") + make_caster<V>::Name + const_name("]"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        if (src.is_none()) {
            value.clear();
            return true;
        }

        if (!PyDict_Check(src.ptr())) {
            return false;
        }

        value.clear();

        PyObject* key;
        PyObject* val;
        Py_ssize_t pos = 0;

        while (PyDict_Next(src.ptr(), &pos, &key, &val)) {
            // Convert key
            OpenMS::String cpp_key;
            if (PyUnicode_Check(key)) {
                Py_ssize_t size;
                const char* data = PyUnicode_AsUTF8AndSize(key, &size);
                if (!data) {
                    PyErr_Clear();
                    return false;
                }
                cpp_key = OpenMS::String(data, size);
            } else {
                return false;
            }

            // Convert value
            type_caster<V> value_caster;
            if (!value_caster.from_python(handle(val), flags, cleanup)) {
                return false;
            }

            value[cpp_key] = std::move(value_caster.value);
        }

        return true;
    }

    static handle from_cpp(const MapType& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        PyObject* dict = PyDict_New();
        if (!dict) return handle();

        for (const auto& kv : src) {
            PyObject* key = PyUnicode_FromStringAndSize(kv.first.c_str(), kv.first.size());
            if (!key) {
                Py_DECREF(dict);
                return handle();
            }

            handle val_handle = type_caster<V>::from_cpp(kv.second, policy, cleanup);
            if (!val_handle) {
                Py_DECREF(key);
                Py_DECREF(dict);
                return handle();
            }

            if (PyDict_SetItem(dict, key, val_handle.ptr()) == -1) {
                Py_DECREF(val_handle.ptr());
                Py_DECREF(key);
                Py_DECREF(dict);
                return handle();
            }

            Py_DECREF(val_handle.ptr());
            Py_DECREF(key);
        }

        return dict;
    }

    static handle from_cpp(MapType& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const MapType&>(src), policy, cleanup);
    }

    static handle from_cpp(MapType&& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const MapType&>(src), policy, cleanup);
    }
};

/**
 * Type caster for OpenMS::StringList (alias for std::vector<OpenMS::String>)
 *
 * Note: StringList is a typedef, so this just delegates to the vector caster.
 * This caster is provided for explicit type matching.
 */
// StringList is already handled by std::vector<OpenMS::String>

/**
 * Type caster for OpenMS::IntList (std::vector<int>)
 */
// IntList uses standard std::vector<int> which nanobind handles natively

/**
 * Type caster for OpenMS::DoubleList (std::vector<double>)
 */
// DoubleList uses standard std::vector<double> which nanobind handles natively

/**
 * Type caster for CVTermMap
 *
 * CVTermMap is: std::map<String, std::vector<CVTerm>>
 * This nested structure needs special handling.
 */
// Note: CVTerm is a class that will be wrapped normally.
// The map caster for std::map<OpenMS::String, V> handles this when V = std::vector<CVTerm>

}  // namespace detail
}  // namespace nanobind

/**
 * Convenience typedefs for common container types
 */
namespace pyopenms {

using StringList = std::vector<OpenMS::String>;
using IntList = std::vector<int>;
using DoubleList = std::vector<double>;

}  // namespace pyopenms
