/**
 * @file limited_api_compat.h
 * @brief Shims for CPython APIs that differ under the Limited API (abi3).
 *
 * Built with PYOPENMS_STABLE_ABI=ON, the extension modules are compiled with
 * Py_LIMITED_API defined and may only use the stable subset of the CPython C
 * API. A few of the fast accessor *macros* are outside that subset and have to
 * fall back to their function counterparts.
 */

#pragma once

#include <Python.h>

/**
 * Store @p item at index @p i of the freshly created list @p list, stealing the
 * reference to @p item.
 *
 * PyList_SET_ITEM() is a macro that writes straight into the list's item array
 * and is not part of the limited API. PyList_SetItem() is, and also steals the
 * reference, but it is an out-of-line call that additionally bounds-checks and
 * releases any previous occupant. That costs about 12% on list-valued DataValue
 * conversions (measured on the "DataValue IntList round-trip" benchmark), so
 * keep the macro wherever it is available.
 *
 * Only use this on a list that was just built with PyList_New(n) and only with
 * i < n, so that the discarded bounds check cannot fail.
 */
#if defined(Py_LIMITED_API)
#  define PYOPENMS_LIST_SET_ITEM(list, i, item) PyList_SetItem((list), (i), (item))
#else
#  define PYOPENMS_LIST_SET_ITEM(list, i, item) PyList_SET_ITEM((list), (i), (item))
#endif
