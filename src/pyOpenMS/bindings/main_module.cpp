// pyOpenMS nanobind bindings
// Main module - placeholder, actual bindings in _pyopenms_<domain> modules

#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(_pyopenms, m) {
    m.doc() = "pyOpenMS: Python bindings for OpenMS (nanobind)\n"
             "Note: Import pyopenms instead for all classes.";
}