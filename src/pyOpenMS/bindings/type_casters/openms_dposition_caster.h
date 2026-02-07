/**
 * @file dposition.h
 * @brief nanobind type casters for OpenMS::DPosition<D>
 *
 * Provides conversion between Python types and OpenMS::DPosition:
 *   - DPosition<1>: float <-> double
 *   - DPosition<2>: tuple/list (x, y) <-> (rt, mz)
 *   - std::vector<DPosition<2>>: numpy 2D array
 */

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <vector>

namespace nanobind {
namespace detail {

/**
 * Type caster for DPosition<1> (1D position)
 *
 * Python -> C++: float/int -> DPosition<1>
 * C++ -> Python: DPosition<1> -> float
 */
template <>
struct type_caster<OpenMS::DPosition<1>> {
public:
    NB_TYPE_CASTER(OpenMS::DPosition<1>, const_name("float"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        if (src.is_none()) {
            return false;
        }

        // Try to convert to double
        PyObject* float_obj = PyNumber_Float(src.ptr());
        if (!float_obj) {
            PyErr_Clear();
            return false;
        }

        double val = PyFloat_AsDouble(float_obj);
        Py_DECREF(float_obj);

        if (PyErr_Occurred()) {
            PyErr_Clear();
            return false;
        }

        value[0] = val;
        return true;
    }

    static handle from_cpp(const OpenMS::DPosition<1>& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return PyFloat_FromDouble(src[0]);
    }

    static handle from_cpp(OpenMS::DPosition<1>& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::DPosition<1>&>(src), policy, cleanup);
    }

    static handle from_cpp(OpenMS::DPosition<1>&& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::DPosition<1>&>(src), policy, cleanup);
    }
};

/**
 * Type caster for DPosition<2> (2D position - RT, MZ)
 *
 * Python -> C++: tuple/list (x, y) -> DPosition<2>
 * C++ -> Python: DPosition<2> -> tuple (rt, mz)
 */
template <>
struct type_caster<OpenMS::DPosition<2>> {
public:
    NB_TYPE_CASTER(OpenMS::DPosition<2>, const_name("tuple[float, float]"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        if (src.is_none()) {
            return false;
        }

        // Must be a sequence of length 2
        if (!PySequence_Check(src.ptr())) {
            return false;
        }

        Py_ssize_t len = PySequence_Length(src.ptr());
        if (len != 2) {
            return false;
        }

        // Extract elements
        PyObject* item0 = PySequence_GetItem(src.ptr(), 0);
        PyObject* item1 = PySequence_GetItem(src.ptr(), 1);

        if (!item0 || !item1) {
            Py_XDECREF(item0);
            Py_XDECREF(item1);
            PyErr_Clear();
            return false;
        }

        // Convert to doubles
        PyObject* float0 = PyNumber_Float(item0);
        PyObject* float1 = PyNumber_Float(item1);

        Py_DECREF(item0);
        Py_DECREF(item1);

        if (!float0 || !float1) {
            Py_XDECREF(float0);
            Py_XDECREF(float1);
            PyErr_Clear();
            return false;
        }

        value[0] = PyFloat_AsDouble(float0);
        value[1] = PyFloat_AsDouble(float1);

        Py_DECREF(float0);
        Py_DECREF(float1);

        if (PyErr_Occurred()) {
            PyErr_Clear();
            return false;
        }

        return true;
    }

    static handle from_cpp(const OpenMS::DPosition<2>& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        PyObject* f0 = PyFloat_FromDouble(src[0]);
        if (!f0) return handle();
        PyObject* f1 = PyFloat_FromDouble(src[1]);
        if (!f1) { Py_DECREF(f0); return handle(); }
        PyObject* tup = PyTuple_Pack(2, f0, f1);
        Py_DECREF(f0);
        Py_DECREF(f1);
        return tup;
    }

    static handle from_cpp(OpenMS::DPosition<2>& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::DPosition<2>&>(src), policy, cleanup);
    }

    static handle from_cpp(OpenMS::DPosition<2>&& src, rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const OpenMS::DPosition<2>&>(src), policy, cleanup);
    }
};

/**
 * Type caster for std::vector<DPosition<2>>
 *
 * Python -> C++: numpy 2D array (Nx2) or list of tuples
 * C++ -> Python: numpy 2D array (Nx2)
 */
template <>
struct type_caster<std::vector<OpenMS::DPosition<2>>> {
public:
    NB_TYPE_CASTER(std::vector<OpenMS::DPosition<2>>,
                   const_name("numpy.ndarray[numpy.float64[Any, 2]]"))

    bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
        if (src.is_none()) {
            value.clear();
            return true;
        }

        // Try numpy array first
        try {
            auto arr = cast<ndarray<double, shape<-1, 2>, c_contig>>(src);
            size_t n = arr.shape(0);
            value.resize(n);

            const double* data = arr.data();
            for (size_t i = 0; i < n; ++i) {
                value[i][0] = data[i * 2];
                value[i][1] = data[i * 2 + 1];
            }
            return true;
        } catch (...) {
            // Fall through to sequence conversion
        }

        // Try as sequence of sequences
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

        type_caster<OpenMS::DPosition<2>> inner_caster;
        for (Py_ssize_t i = 0; i < len; ++i) {
            PyObject* item = PySequence_GetItem(src.ptr(), i);
            if (!item) {
                PyErr_Clear();
                return false;
            }

            handle item_handle(item);
            if (!inner_caster.from_python(item_handle, flags, cleanup)) {
                Py_DECREF(item);
                return false;
            }

            value.push_back(inner_caster.value);
            Py_DECREF(item);
        }

        return true;
    }

    static handle from_cpp(const std::vector<OpenMS::DPosition<2>>& src,
                          rv_policy policy, cleanup_list* cleanup) noexcept {
        size_t n = src.size();

        // Create Python list of tuples
        PyObject* list = PyList_New(n);
        if (!list) return handle();

        for (size_t i = 0; i < n; ++i) {
            PyObject* f0 = PyFloat_FromDouble(src[i][0]);
            if (!f0) { Py_DECREF(list); return handle(); }
            PyObject* f1 = PyFloat_FromDouble(src[i][1]);
            if (!f1) { Py_DECREF(f0); Py_DECREF(list); return handle(); }
            PyObject* tuple = PyTuple_Pack(2, f0, f1);
            Py_DECREF(f0);
            Py_DECREF(f1);
            if (!tuple) {
                Py_DECREF(list);
                return handle();
            }
            PyList_SET_ITEM(list, i, tuple);  // steals reference
        }

        return handle(list);
    }

    static handle from_cpp(std::vector<OpenMS::DPosition<2>>& src,
                          rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const std::vector<OpenMS::DPosition<2>>&>(src), policy, cleanup);
    }

    static handle from_cpp(std::vector<OpenMS::DPosition<2>>&& src,
                          rv_policy policy, cleanup_list* cleanup) noexcept {
        return from_cpp(const_cast<const std::vector<OpenMS::DPosition<2>>&>(src), policy, cleanup);
    }
};

// Type aliases for convenience
using DPosition1 = OpenMS::DPosition<1>;
using DPosition2 = OpenMS::DPosition<2>;

}  // namespace detail
}  // namespace nanobind
