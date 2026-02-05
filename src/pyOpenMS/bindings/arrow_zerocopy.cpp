// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2024.
//
// This software is released under a three-clause BSD license.
// --------------------------------------------------------------------------
// Zero-copy Arrow export for pyOpenMS via nanobind
// --------------------------------------------------------------------------

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <arrow/c/abi.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>

#include <cstdlib>
#include <cstring>

namespace nb = nanobind;

// Helper: RAII guard for malloc'd Arrow structs
struct ArrowGuard {
    ArrowSchema* schema;
    ArrowArray* array;

    ArrowGuard()
        : schema(static_cast<ArrowSchema*>(std::malloc(sizeof(ArrowSchema)))),
          array(static_cast<ArrowArray*>(std::malloc(sizeof(ArrowArray))))
    {
        if (!schema || !array) {
            std::free(schema);
            std::free(array);
            throw std::bad_alloc();
        }
        std::memset(schema, 0, sizeof(ArrowSchema));
        std::memset(array, 0, sizeof(ArrowArray));
    }

    ~ArrowGuard() {
        if (schema) {
            if (schema->release) schema->release(schema);
            std::free(schema);
        }
        if (array) {
            if (array->release) array->release(array);
            std::free(array);
        }
    }

    // Release ownership (after PyArrow takes over)
    void release() {
        // PyArrow now owns the Arrow data; just free the malloc'd structs
        // without calling the release callbacks
        std::free(schema);
        std::free(array);
        schema = nullptr;
        array = nullptr;
    }

    ArrowGuard(const ArrowGuard&) = delete;
    ArrowGuard& operator=(const ArrowGuard&) = delete;
};

// Import a filled ArrowSchema+ArrowArray into a PyArrow Table
static nb::object import_to_pyarrow(ArrowGuard& guard) {
    nb::module_ pa = nb::module_::import_("pyarrow");

    // _import_from_c expects integer addresses
    auto batch = pa.attr("RecordBatch").attr("_import_from_c")(
        reinterpret_cast<uintptr_t>(guard.array),
        reinterpret_cast<uintptr_t>(guard.schema)
    );

    // PyArrow now owns the Arrow data
    guard.release();

    return pa.attr("Table").attr("from_batches")(nb::make_tuple(batch));
}

NB_MODULE(_arrow_zerocopy, m) {
    m.doc() = "Zero-copy Arrow export from OpenMS MSExperiment.\n\n"
              "This module provides zero-copy export of MS data to Apache Arrow format "
              "using the Arrow C Data Interface. It is only available when OpenMS is "
              "built with WITH_PARQUET=ON.\n\n"
              ".. warning::\n"
              "    **EXPERIMENTAL API**: This module is experimental and may change.";

    m.def("spectra_to_arrow",
        [](nb::object exp_obj,
           std::string format,
           nb::object ms_levels_obj,
           double min_rt,
           double max_rt,
           double min_mz,
           double max_mz,
           nb::object columns_obj,
           bool include_precursor_info,
           bool include_ion_mobility) -> nb::object
        {
            // Cast to C++ reference via nanobind
            const auto& exp = nb::cast<const OpenMS::MSExperiment&>(exp_obj);

            // Build config
            OpenMS::ArrowSpectraExportConfig config;

            // Format
            std::string fmt = format;
            for (auto& c : fmt) c = static_cast<char>(std::tolower(c));
            if (fmt == "long")
                config.format = OpenMS::ArrowExportFormat::Long;
            else if (fmt == "semi_wide")
                config.format = OpenMS::ArrowExportFormat::SemiWide;
            else
                throw nb::value_error(("format must be 'long' or 'semi_wide', got '" + format + "'").c_str());

            config.min_rt = min_rt;
            config.max_rt = max_rt;
            config.min_mz = min_mz;
            config.max_mz = max_mz;
            config.include_precursor_info = include_precursor_info;
            config.include_ion_mobility = include_ion_mobility;

            // ms_levels
            if (!ms_levels_obj.is_none()) {
                for (auto item : ms_levels_obj) {
                    config.ms_levels.push_back(nb::cast<unsigned int>(item));
                }
            }

            // columns
            if (!columns_obj.is_none()) {
                for (auto item : columns_obj) {
                    config.columns.push_back(nb::cast<std::string>(item));
                }
            }

            // Allocate Arrow structs
            ArrowGuard guard;

            // Call C++ export (release GIL)
            bool success;
            {
                nb::gil_scoped_release release;
                success = OpenMS::MSExperimentArrowExport::exportSpectraToArrowCDataInterface(
                    exp, config, guard.schema, guard.array);
            }

            if (!success)
                throw std::runtime_error("Failed to export spectra to Arrow format");

            return import_to_pyarrow(guard);
        },
        nb::arg("exp"),
        nb::arg("format") = "long",
        nb::arg("ms_levels") = nb::none(),
        nb::arg("min_rt") = 0.0,
        nb::arg("max_rt") = 0.0,
        nb::arg("min_mz") = 0.0,
        nb::arg("max_mz") = 0.0,
        nb::arg("columns") = nb::none(),
        nb::arg("include_precursor_info") = true,
        nb::arg("include_ion_mobility") = true,
        "Export spectra to Arrow Table using zero-copy C Data Interface.\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("chromatograms_to_arrow",
        [](nb::object exp_obj,
           std::string format,
           double min_rt,
           double max_rt,
           nb::object columns_obj) -> nb::object
        {
            const auto& exp = nb::cast<const OpenMS::MSExperiment&>(exp_obj);

            OpenMS::ArrowChromatogramExportConfig config;

            std::string fmt = format;
            for (auto& c : fmt) c = static_cast<char>(std::tolower(c));
            if (fmt == "long")
                config.format = OpenMS::ArrowExportFormat::Long;
            else if (fmt == "semi_wide")
                config.format = OpenMS::ArrowExportFormat::SemiWide;
            else
                throw nb::value_error(("format must be 'long' or 'semi_wide', got '" + format + "'").c_str());

            config.min_rt = min_rt;
            config.max_rt = max_rt;

            if (!columns_obj.is_none()) {
                for (auto item : columns_obj) {
                    config.columns.push_back(nb::cast<std::string>(item));
                }
            }

            ArrowGuard guard;

            bool success;
            {
                nb::gil_scoped_release release;
                success = OpenMS::MSExperimentArrowExport::exportChromatogramsToArrowCDataInterface(
                    exp, config, guard.schema, guard.array);
            }

            if (!success)
                throw std::runtime_error("Failed to export chromatograms to Arrow format");

            return import_to_pyarrow(guard);
        },
        nb::arg("exp"),
        nb::arg("format") = "long",
        nb::arg("min_rt") = 0.0,
        nb::arg("max_rt") = 0.0,
        nb::arg("columns") = nb::none(),
        "Export chromatograms to Arrow Table using zero-copy C Data Interface.\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );
}
