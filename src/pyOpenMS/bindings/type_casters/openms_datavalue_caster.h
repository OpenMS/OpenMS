/**
 * @file datavalue.h
 * @brief nanobind type casters for OpenMS::DataValue and OpenMS::ParamValue
 *
 * These are variant types that can hold multiple value types.
 * DataValue: int, double, String, StringList, IntList, DoubleList
 * ParamValue: similar to DataValue with slight differences
 */

#pragma once

#include <nanobind/nanobind.h>
// Use custom std::string caster that accepts bytes (do NOT include nanobind/stl/string.h)
#include "std_string_bytes_caster.h"
#include <nanobind/stl/vector.h>
#include <limits>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include "string.h"  // For OpenMS::String caster

namespace nanobind {
namespace detail {

/**
 * Type caster for OpenMS::DataValue
 *
 * DataValue is a variant type that can hold:
 *   - Empty value
 *   - int (Int64)
 *   - double
 *   - String
 *   - StringList (vector<String>)
 *   - IntList (vector<int>)
 *   - DoubleList (vector<double>)
 *
 * Python -> C++:
 *   - None -> Empty
 *   - int -> Int
 *   - float -> Double
 *   - str -> String
 *   - bytes -> String
 *   - list of str -> StringList
 *   - list of int -> IntList
 *   - list of float -> DoubleList
 *
 * C++ -> Python:
 *   - Converts to appropriate Python type
 */
template <>
struct type_caster<OpenMS::DataValue> {
public:
    NB_TYPE_CASTER(OpenMS::DataValue,
                   const_name("None | int | float | str | bytes | "
                             "list[str] | list[int] | list[float]"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        // Handle None -> Empty DataValue
        if (src.is_none()) {
            value = OpenMS::DataValue();
            return true;
        }

        // Try int (but not bool, which is a subclass of int in Python)
        if (PyLong_Check(src.ptr()) && !PyBool_Check(src.ptr())) {
            long long val = PyLong_AsLongLong(src.ptr());
            if (PyErr_Occurred()) {
                PyErr_Clear();
                return false;
            }
            value = OpenMS::DataValue(static_cast<int64_t>(val));
            return true;
        }

        // Try float
        if (PyFloat_Check(src.ptr())) {
            double val = PyFloat_AsDouble(src.ptr());
            if (PyErr_Occurred()) {
                PyErr_Clear();
                return false;
            }
            value = OpenMS::DataValue(val);
            return true;
        }

        // Try str
        if (PyUnicode_Check(src.ptr())) {
            Py_ssize_t size;
            const char* data = PyUnicode_AsUTF8AndSize(src.ptr(), &size);
            if (!data) {
                PyErr_Clear();
                return false;
            }
            value = OpenMS::DataValue(OpenMS::String(data, size));
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
            value = OpenMS::DataValue(OpenMS::String(data, size));
            return true;
        }

        // Try list
        if (PyList_Check(src.ptr()) || PyTuple_Check(src.ptr())) {
            Py_ssize_t len = PySequence_Length(src.ptr());
            if (len < 0) {
                PyErr_Clear();
                return false;
            }

            if (len == 0) {
                // Empty list - default to StringList
                value = OpenMS::DataValue(OpenMS::StringList());
                return true;
            }

            // Check first element to determine type
            PyObject* first = PySequence_GetItem(src.ptr(), 0);
            if (!first) {
                PyErr_Clear();
                return false;
            }

            bool is_string = PyUnicode_Check(first) || PyBytes_Check(first);
            bool is_int = PyLong_Check(first) && !PyBool_Check(first);
            bool is_float = PyFloat_Check(first);

            Py_DECREF(first);

            if (is_string) {
                // StringList
                OpenMS::StringList sl;
                sl.reserve(len);
                for (Py_ssize_t i = 0; i < len; ++i) {
                    PyObject* item = PySequence_GetItem(src.ptr(), i);
                    if (!item) {
                        PyErr_Clear();
                        return false;
                    }

                    if (PyUnicode_Check(item)) {
                        Py_ssize_t size;
                        const char* data = PyUnicode_AsUTF8AndSize(item, &size);
                        if (!data) {
                            Py_DECREF(item);
                            PyErr_Clear();
                            return false;
                        }
                        sl.push_back(OpenMS::String(data, size));
                    } else if (PyBytes_Check(item)) {
                        char* data;
                        Py_ssize_t size;
                        if (PyBytes_AsStringAndSize(item, &data, &size) != 0) {
                            Py_DECREF(item);
                            PyErr_Clear();
                            return false;
                        }
                        sl.push_back(OpenMS::String(data, size));
                    } else {
                        Py_DECREF(item);
                        return false;  // reject non-string items
                    }
                    Py_DECREF(item);
                }
                value = OpenMS::DataValue(sl);
                return true;
            } else if (is_int) {
                // IntList
                OpenMS::IntList il;
                il.reserve(len);
                for (Py_ssize_t i = 0; i < len; ++i) {
                    PyObject* item = PySequence_GetItem(src.ptr(), i);
                    if (!item) {
                        PyErr_Clear();
                        return false;
                    }
                    long val = PyLong_AsLong(item);
                    Py_DECREF(item);
                    if (PyErr_Occurred()) {
                        PyErr_Clear();
                        return false;
                    }
                    if (val > std::numeric_limits<int>::max() || val < std::numeric_limits<int>::min()) {
                        return false;  // overflow
                    }
                    il.push_back(static_cast<int>(val));
                }
                value = OpenMS::DataValue(il);
                return true;
            } else if (is_float) {
                // DoubleList
                OpenMS::DoubleList dl;
                dl.reserve(len);
                for (Py_ssize_t i = 0; i < len; ++i) {
                    PyObject* item = PySequence_GetItem(src.ptr(), i);
                    if (!item) {
                        PyErr_Clear();
                        return false;
                    }
                    PyObject* float_item = PyNumber_Float(item);
                    Py_DECREF(item);
                    if (!float_item) {
                        PyErr_Clear();
                        return false;
                    }
                    double val = PyFloat_AsDouble(float_item);
                    Py_DECREF(float_item);
                    dl.push_back(val);
                }
                value = OpenMS::DataValue(dl);
                return true;
            }
        }

        return false;
    }

    static handle from_cpp(const OpenMS::DataValue& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        using ValueType = OpenMS::DataValue::DataType;

        switch (src.valueType()) {
            case ValueType::EMPTY_VALUE:
                return none().release();

            case ValueType::INT_VALUE:
                return PyLong_FromLongLong(static_cast<long long>(src));

            case ValueType::DOUBLE_VALUE:
                return PyFloat_FromDouble(static_cast<double>(src));

            case ValueType::STRING_VALUE: {
                // Use toString() - DataValue has no operator String()
                OpenMS::String s = src.toString();
                return PyUnicode_FromStringAndSize(s.c_str(), s.size());
            }

            case ValueType::STRING_LIST: {
                const OpenMS::StringList& sl = static_cast<const OpenMS::StringList&>(src);
                PyObject* list = PyList_New(sl.size());
                if (!list) return handle();
                for (size_t i = 0; i < sl.size(); ++i) {
                    PyObject* item = PyUnicode_FromStringAndSize(sl[i].c_str(), sl[i].size());
                    if (!item) {
                        Py_DECREF(list);
                        return handle();
                    }
                    PyList_SET_ITEM(list, i, item);
                }
                return list;
            }

            case ValueType::INT_LIST: {
                const OpenMS::IntList& il = static_cast<const OpenMS::IntList&>(src);
                PyObject* list = PyList_New(il.size());
                if (!list) return handle();
                for (size_t i = 0; i < il.size(); ++i) {
                    PyObject* item = PyLong_FromLong(il[i]);
                    if (!item) {
                        Py_DECREF(list);
                        return handle();
                    }
                    PyList_SET_ITEM(list, i, item);
                }
                return list;
            }

            case ValueType::DOUBLE_LIST: {
                const OpenMS::DoubleList& dl = static_cast<const OpenMS::DoubleList&>(src);
                PyObject* list = PyList_New(dl.size());
                if (!list) return handle();
                for (size_t i = 0; i < dl.size(); ++i) {
                    PyObject* item = PyFloat_FromDouble(dl[i]);
                    if (!item) {
                        Py_DECREF(list);
                        return handle();
                    }
                    PyList_SET_ITEM(list, i, item);
                }
                return list;
            }

            default:
                return none().release();
        }
    }

    static handle from_cpp(OpenMS::DataValue& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::DataValue&>(src), policy, cleanup);
    }

    static handle from_cpp(OpenMS::DataValue&& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::DataValue&>(src), policy, cleanup);
    }
};

/**
 * Type caster for OpenMS::ParamValue
 *
 * ParamValue is similar to DataValue but used in Param objects.
 * Conversion logic is essentially the same.
 */
template <>
struct type_caster<OpenMS::ParamValue> {
public:
    NB_TYPE_CASTER(OpenMS::ParamValue,
                   const_name("None | int | float | str | bytes | "
                             "list[str] | list[int] | list[float]"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        // Handle None -> Empty ParamValue
        if (src.is_none()) {
            value = OpenMS::ParamValue();
            return true;
        }

        // Try int (but not bool)
        if (PyLong_Check(src.ptr()) && !PyBool_Check(src.ptr())) {
            long long val = PyLong_AsLongLong(src.ptr());
            if (PyErr_Occurred()) {
                PyErr_Clear();
                return false;
            }
            value = OpenMS::ParamValue(static_cast<int64_t>(val));
            return true;
        }

        // Try float
        if (PyFloat_Check(src.ptr())) {
            double val = PyFloat_AsDouble(src.ptr());
            if (PyErr_Occurred()) {
                PyErr_Clear();
                return false;
            }
            value = OpenMS::ParamValue(val);
            return true;
        }

        // Try str
        if (PyUnicode_Check(src.ptr())) {
            Py_ssize_t size;
            const char* data = PyUnicode_AsUTF8AndSize(src.ptr(), &size);
            if (!data) {
                PyErr_Clear();
                return false;
            }
            value = OpenMS::ParamValue(std::string(data, size));
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
            value = OpenMS::ParamValue(std::string(data, size));
            return true;
        }

        // Try list - same logic as DataValue
        if (PyList_Check(src.ptr()) || PyTuple_Check(src.ptr())) {
            Py_ssize_t len = PySequence_Length(src.ptr());
            if (len < 0) {
                PyErr_Clear();
                return false;
            }

            if (len == 0) {
                value = OpenMS::ParamValue(std::vector<std::string>());
                return true;
            }

            PyObject* first = PySequence_GetItem(src.ptr(), 0);
            if (!first) {
                PyErr_Clear();
                return false;
            }

            bool is_string = PyUnicode_Check(first) || PyBytes_Check(first);
            bool is_int = PyLong_Check(first) && !PyBool_Check(first);
            bool is_float = PyFloat_Check(first);

            Py_DECREF(first);

            if (is_string) {
                std::vector<std::string> sl;
                sl.reserve(len);
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
                            sl.push_back(std::string(data, size));
                        }
                    } else if (PyBytes_Check(item)) {
                        char* data;
                        Py_ssize_t size;
                        if (PyBytes_AsStringAndSize(item, &data, &size) == 0) {
                            sl.push_back(std::string(data, size));
                        }
                    }
                    Py_DECREF(item);
                }
                value = OpenMS::ParamValue(sl);
                return true;
            } else if (is_int) {
                std::vector<int> il;
                il.reserve(len);
                for (Py_ssize_t i = 0; i < len; ++i) {
                    PyObject* item = PySequence_GetItem(src.ptr(), i);
                    if (!item) {
                        PyErr_Clear();
                        return false;
                    }
                    long val = PyLong_AsLong(item);
                    Py_DECREF(item);
                    if (PyErr_Occurred()) {
                        PyErr_Clear();
                        return false;
                    }
                    if (val > std::numeric_limits<int>::max() || val < std::numeric_limits<int>::min()) {
                        return false;  // overflow
                    }
                    il.push_back(static_cast<int>(val));
                }
                value = OpenMS::ParamValue(il);
                return true;
            } else if (is_float) {
                std::vector<double> dl;
                dl.reserve(len);
                for (Py_ssize_t i = 0; i < len; ++i) {
                    PyObject* item = PySequence_GetItem(src.ptr(), i);
                    if (!item) {
                        PyErr_Clear();
                        return false;
                    }
                    PyObject* float_item = PyNumber_Float(item);
                    Py_DECREF(item);
                    if (!float_item) {
                        PyErr_Clear();
                        return false;
                    }
                    double val = PyFloat_AsDouble(float_item);
                    Py_DECREF(float_item);
                    dl.push_back(val);
                }
                value = OpenMS::ParamValue(dl);
                return true;
            }
        }

        return false;
    }

    static handle from_cpp(const OpenMS::ParamValue& src, rv_policy policy,
                          cleanup_list* cleanup) noexcept {
        using ValueType = OpenMS::ParamValue::ValueType;

        switch (src.valueType()) {
            case ValueType::EMPTY_VALUE:
                return none().release();

            case ValueType::INT_VALUE:
                return PyLong_FromLongLong(static_cast<long long>(static_cast<int>(src)));  // ParamValue stores int

            case ValueType::DOUBLE_VALUE:
                return PyFloat_FromDouble(static_cast<double>(src));

            case ValueType::STRING_VALUE: {
                const std::string s = src;
                return PyUnicode_FromStringAndSize(s.c_str(), s.size());
            }

            case ValueType::STRING_LIST: {
                const std::vector<std::string> sl = src;
                PyObject* list = PyList_New(sl.size());
                if (!list) return handle();
                for (size_t i = 0; i < sl.size(); ++i) {
                    PyObject* item = PyUnicode_FromStringAndSize(sl[i].c_str(), sl[i].size());
                    if (!item) {
                        Py_DECREF(list);
                        return handle();
                    }
                    PyList_SET_ITEM(list, i, item);
                }
                return list;
            }

            case ValueType::INT_LIST: {
                const std::vector<int> il = src;
                PyObject* list = PyList_New(il.size());
                if (!list) return handle();
                for (size_t i = 0; i < il.size(); ++i) {
                    PyObject* item = PyLong_FromLong(il[i]);
                    if (!item) {
                        Py_DECREF(list);
                        return handle();
                    }
                    PyList_SET_ITEM(list, i, item);
                }
                return list;
            }

            case ValueType::DOUBLE_LIST: {
                const std::vector<double> dl = src;
                PyObject* list = PyList_New(dl.size());
                if (!list) return handle();
                for (size_t i = 0; i < dl.size(); ++i) {
                    PyObject* item = PyFloat_FromDouble(dl[i]);
                    if (!item) {
                        Py_DECREF(list);
                        return handle();
                    }
                    PyList_SET_ITEM(list, i, item);
                }
                return list;
            }

            default:
                return none().release();
        }
    }

    static handle from_cpp(OpenMS::ParamValue& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::ParamValue&>(src), policy, cleanup);
    }

    static handle from_cpp(OpenMS::ParamValue&& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::ParamValue&>(src), policy, cleanup);
    }
};

}  // namespace detail
}  // namespace nanobind
