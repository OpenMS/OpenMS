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
#include <arrow/c/bridge.h>
#include <arrow/api.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>
#include <OpenMS/FORMAT/FeatureMapArrowIO.h>
#include <OpenMS/FORMAT/ConsensusMapArrowIO.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>

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

// Convert a C++ arrow::Table to a PyArrow Table via the C Data Interface.
// The table is combined into a single RecordBatch, exported to C structs,
// then imported by PyArrow (zero-copy).
static nb::object table_to_pyarrow(const std::shared_ptr<arrow::Table>& table) {
    if (!table)
        throw std::runtime_error("Arrow export returned null table");

    auto batch_result = table->CombineChunksToBatch();
    if (!batch_result.ok())
        throw std::runtime_error("Failed to combine Arrow table chunks: " + batch_result.status().ToString());

    auto batch = batch_result.ValueOrDie();

    ArrowGuard guard;
    auto schema_status = arrow::ExportSchema(*batch->schema(), guard.schema);
    if (!schema_status.ok())
        throw std::runtime_error("Failed to export Arrow schema: " + schema_status.ToString());

    auto array_status = arrow::ExportRecordBatch(*batch, guard.array);
    if (!array_status.ok())
        throw std::runtime_error("Failed to export Arrow record batch: " + array_status.ToString());

    return import_to_pyarrow(guard);
}

// Convert a PyArrow Table to a C++ arrow::Table via the C Data Interface.
static std::shared_ptr<arrow::Table> pyarrow_to_table(nb::object pa_table) {
    nb::module_ pa = nb::module_::import_("pyarrow");

    // Ensure we have a Table
    if (!nb::isinstance(pa_table, pa.attr("Table")))
        throw nb::type_error("Expected a pyarrow.Table");

    // Combine chunks into a single RecordBatch for C Data Interface export
    nb::object combined = pa_table.attr("combine_chunks")();
    auto batches = nb::cast<nb::list>(combined.attr("to_batches")());
    if (batches.size() == 0)
        throw std::runtime_error("PyArrow table has no record batches");
    nb::object batch = batches[0];

    // Allocate C structs for export
    ArrowGuard guard;

    // Export from PyArrow to C Data Interface
    batch.attr("_export_to_c")(
        reinterpret_cast<uintptr_t>(guard.array),
        reinterpret_cast<uintptr_t>(guard.schema)
    );

    // Import into C++ Arrow (consumes the C structs — release callbacks become null)
    auto schema_result = arrow::ImportSchema(guard.schema);
    if (!schema_result.ok())
        throw std::runtime_error("Failed to import Arrow schema: " + schema_result.status().ToString());

    auto batch_result = arrow::ImportRecordBatch(guard.array, *schema_result);
    if (!batch_result.ok())
        throw std::runtime_error("Failed to import Arrow record batch: " + batch_result.status().ToString());

    // C structs have been consumed by Import (release callbacks are now null),
    // so the guard destructor will just free the malloc'd memory
    return arrow::Table::Make((*batch_result)->schema(),
                              (*batch_result)->columns(),
                              (*batch_result)->num_rows());
}

NB_MODULE(_arrow_zerocopy, m) {
    m.doc() = "Zero-copy Arrow export from OpenMS MSExperiment.\n\n"
              "This module provides zero-copy export of MS data to Apache Arrow format "
              "using the Arrow C Data Interface.\n\n"
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

    // -----------------------------------------------------------------------
    // FeatureMapArrowIO — zero-copy Arrow export/import
    // -----------------------------------------------------------------------

    m.def("featuremap_features_to_arrow",
        [](nb::object fmap_obj) -> nb::object
        {
            const auto& fmap = nb::cast<const OpenMS::FeatureMap&>(fmap_obj);
            std::shared_ptr<arrow::Table> table;
            {
                nb::gil_scoped_release release;
                table = OpenMS::FeatureMapArrowIO::exportFeaturesToArrow(fmap);
            }
            return table_to_pyarrow(table);
        },
        nb::arg("feature_map"),
        "Export FeatureMap features to a PyArrow Table (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("featuremap_psms_to_arrow",
        [](nb::object fmap_obj) -> nb::object
        {
            const auto& fmap = nb::cast<const OpenMS::FeatureMap&>(fmap_obj);
            std::shared_ptr<arrow::Table> table;
            {
                nb::gil_scoped_release release;
                table = OpenMS::FeatureMapArrowIO::exportPSMsToArrow(fmap);
            }
            return table_to_pyarrow(table);
        },
        nb::arg("feature_map"),
        "Export FeatureMap PSMs to a PyArrow Table (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("featuremap_import_features_from_arrow",
        [](nb::object pa_table, nb::object fmap_obj) -> bool
        {
            auto& fmap = nb::cast<OpenMS::FeatureMap&>(fmap_obj);
            auto table = pyarrow_to_table(pa_table);
            nb::gil_scoped_release release;
            return OpenMS::FeatureMapArrowIO::importFeaturesFromArrow(table, fmap);
        },
        nb::arg("table"), nb::arg("feature_map"),
        "Import features from a PyArrow Table into a FeatureMap (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("featuremap_import_psms_from_arrow",
        [](nb::object pa_table, nb::object fmap_obj) -> bool
        {
            auto& fmap = nb::cast<OpenMS::FeatureMap&>(fmap_obj);
            auto table = pyarrow_to_table(pa_table);
            nb::gil_scoped_release release;
            return OpenMS::FeatureMapArrowIO::importPSMsFromArrow(table, fmap);
        },
        nb::arg("table"), nb::arg("feature_map"),
        "Import PSMs from a PyArrow Table into a FeatureMap (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    // -----------------------------------------------------------------------
    // ConsensusMapArrowIO — zero-copy Arrow export/import
    // -----------------------------------------------------------------------

    m.def("consensusmap_features_to_arrow",
        [](nb::object cmap_obj) -> nb::object
        {
            const auto& cmap = nb::cast<const OpenMS::ConsensusMap&>(cmap_obj);
            std::shared_ptr<arrow::Table> table;
            {
                nb::gil_scoped_release release;
                table = OpenMS::ConsensusMapArrowIO::exportFeaturesToArrow(cmap);
            }
            return table_to_pyarrow(table);
        },
        nb::arg("consensus_map"),
        "Export ConsensusMap features to a PyArrow Table (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("consensusmap_psms_to_arrow",
        [](nb::object cmap_obj) -> nb::object
        {
            const auto& cmap = nb::cast<const OpenMS::ConsensusMap&>(cmap_obj);
            std::shared_ptr<arrow::Table> table;
            {
                nb::gil_scoped_release release;
                table = OpenMS::ConsensusMapArrowIO::exportPSMsToArrow(cmap);
            }
            return table_to_pyarrow(table);
        },
        nb::arg("consensus_map"),
        "Export ConsensusMap PSMs to a PyArrow Table (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("consensusmap_import_features_from_arrow",
        [](nb::object pa_table, nb::object cmap_obj) -> bool
        {
            auto& cmap = nb::cast<OpenMS::ConsensusMap&>(cmap_obj);
            auto table = pyarrow_to_table(pa_table);
            nb::gil_scoped_release release;
            return OpenMS::ConsensusMapArrowIO::importFeaturesFromArrow(table, cmap);
        },
        nb::arg("table"), nb::arg("consensus_map"),
        "Import features from a PyArrow Table into a ConsensusMap (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("consensusmap_import_psms_from_arrow",
        [](nb::object pa_table, nb::object cmap_obj) -> bool
        {
            auto& cmap = nb::cast<OpenMS::ConsensusMap&>(cmap_obj);
            auto table = pyarrow_to_table(pa_table);
            nb::gil_scoped_release release;
            return OpenMS::ConsensusMapArrowIO::importPSMsFromArrow(table, cmap);
        },
        nb::arg("table"), nb::arg("consensus_map"),
        "Import PSMs from a PyArrow Table into a ConsensusMap (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    // -----------------------------------------------------------------------
    // ProteinIdentificationArrowIO — zero-copy Arrow export/import
    // -----------------------------------------------------------------------

    m.def("protein_ids_proteins_to_arrow",
        [](nb::object prot_ids_obj) -> nb::object
        {
            auto prot_ids = nb::cast<std::vector<OpenMS::ProteinIdentification>>(prot_ids_obj);
            std::shared_ptr<arrow::Table> table;
            {
                nb::gil_scoped_release release;
                table = OpenMS::ProteinIdentificationArrowIO::exportProteinsToArrow(prot_ids);
            }
            return table_to_pyarrow(table);
        },
        nb::arg("protein_identifications"),
        "Export protein hits to a PyArrow Table (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("protein_ids_groups_to_arrow",
        [](nb::object prot_ids_obj) -> nb::object
        {
            auto prot_ids = nb::cast<std::vector<OpenMS::ProteinIdentification>>(prot_ids_obj);
            std::shared_ptr<arrow::Table> table;
            {
                nb::gil_scoped_release release;
                table = OpenMS::ProteinIdentificationArrowIO::exportProteinGroupsToArrow(prot_ids);
            }
            return table_to_pyarrow(table);
        },
        nb::arg("protein_identifications"),
        "Export protein groups to a PyArrow Table (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("protein_ids_search_params_to_arrow",
        [](nb::object prot_ids_obj) -> nb::object
        {
            auto prot_ids = nb::cast<std::vector<OpenMS::ProteinIdentification>>(prot_ids_obj);
            std::shared_ptr<arrow::Table> table;
            {
                nb::gil_scoped_release release;
                table = OpenMS::ProteinIdentificationArrowIO::exportSearchParamsToArrow(prot_ids);
            }
            return table_to_pyarrow(table);
        },
        nb::arg("protein_identifications"),
        "Export search parameters to a PyArrow Table (zero-copy).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    // Import methods: std::vector<ProteinIdentification>& is an output param.
    // Since vectors go through nanobind's STL type caster (creates copies),
    // we take by value, modify, and return the result.

    m.def("protein_ids_import_search_params_from_arrow",
        [](nb::object pa_table) -> std::vector<OpenMS::ProteinIdentification>
        {
            auto table = pyarrow_to_table(pa_table);
            std::vector<OpenMS::ProteinIdentification> prot_ids;
            bool ok;
            {
                nb::gil_scoped_release release;
                ok = OpenMS::ProteinIdentificationArrowIO::importSearchParamsFromArrow(table, prot_ids);
            }
            if (!ok) throw std::runtime_error("Failed to import search parameters from Arrow table");
            return prot_ids;
        },
        nb::arg("table"),
        "Import search parameters from a PyArrow Table. Returns list of ProteinIdentifications.\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("protein_ids_import_proteins_from_arrow",
        [](nb::object pa_table, std::vector<OpenMS::ProteinIdentification> prot_ids)
            -> std::vector<OpenMS::ProteinIdentification>
        {
            auto table = pyarrow_to_table(pa_table);
            bool ok;
            {
                nb::gil_scoped_release release;
                ok = OpenMS::ProteinIdentificationArrowIO::importProteinsFromArrow(table, prot_ids);
            }
            if (!ok) throw std::runtime_error("Failed to import proteins from Arrow table");
            return prot_ids;
        },
        nb::arg("table"), nb::arg("protein_identifications"),
        "Import protein hits from a PyArrow Table into existing ProteinIdentifications. Returns updated list.\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );

    m.def("protein_ids_import_groups_from_arrow",
        [](nb::object pa_table, std::vector<OpenMS::ProteinIdentification> prot_ids)
            -> std::vector<OpenMS::ProteinIdentification>
        {
            auto table = pyarrow_to_table(pa_table);
            bool ok;
            {
                nb::gil_scoped_release release;
                ok = OpenMS::ProteinIdentificationArrowIO::importProteinGroupsFromArrow(table, prot_ids);
            }
            if (!ok) throw std::runtime_error("Failed to import protein groups from Arrow table");
            return prot_ids;
        },
        nb::arg("table"), nb::arg("protein_identifications"),
        "Import protein groups from a PyArrow Table into existing ProteinIdentifications. Returns updated list.\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );
}
