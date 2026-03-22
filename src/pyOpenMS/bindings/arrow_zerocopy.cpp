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
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>

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

// Export a C++ arrow::Schema to a pyarrow.Schema via the C Data Interface.
static nb::object schema_to_pyarrow(const std::shared_ptr<arrow::Schema>& schema) {
    // Allocate a C ArrowSchema struct for the C Data Interface
    ArrowSchema* c_schema = static_cast<ArrowSchema*>(std::malloc(sizeof(ArrowSchema)));
    if (!c_schema) throw std::bad_alloc();
    std::memset(c_schema, 0, sizeof(ArrowSchema));

    auto status = arrow::ExportSchema(*schema, c_schema);
    if (!status.ok()) {
        if (c_schema->release) c_schema->release(c_schema);
        std::free(c_schema);
        throw std::runtime_error("Failed to export Arrow schema: " + status.ToString());
    }

    nb::module_ pa = nb::module_::import_("pyarrow");
    // PyArrow's _import_from_c consumes the C struct (sets release to null)
    nb::object result = pa.attr("Schema").attr("_import_from_c")(
        reinterpret_cast<uintptr_t>(c_schema));

    // PyArrow has consumed the struct; just free the malloc'd memory
    std::free(c_schema);
    return result;
}

// Import a pyarrow.Schema to a C++ arrow::Schema via the C Data Interface.
static std::shared_ptr<arrow::Schema> pyarrow_schema_to_arrow(nb::object py_schema) {
    // Allocate a C ArrowSchema struct for the C Data Interface
    ArrowSchema* c_schema = static_cast<ArrowSchema*>(std::malloc(sizeof(ArrowSchema)));
    if (!c_schema) throw std::bad_alloc();
    std::memset(c_schema, 0, sizeof(ArrowSchema));

    // Export from PyArrow to C struct
    py_schema.attr("_export_to_c")(reinterpret_cast<uintptr_t>(c_schema));

    // Import into C++ Arrow
    auto result = arrow::ImportSchema(c_schema);

    // C struct has been consumed by ImportSchema; free the malloc'd memory
    std::free(c_schema);

    if (!result.ok())
        throw std::runtime_error("Failed to import Arrow schema: " + result.status().ToString());

    return result.MoveValueUnsafe();
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

    // -----------------------------------------------------------------------
    // ArrowSchemaRegistry — schema struct bindings
    // -----------------------------------------------------------------------

    nb::class_<OpenMS::ProteinSchema>(m, "ProteinSchema")
        .def_prop_ro_static("ACCESSION", [](nb::handle) { return OpenMS::ProteinSchema::ACCESSION; })
        .def_prop_ro_static("SCORE", [](nb::handle) { return OpenMS::ProteinSchema::SCORE; })
        .def_prop_ro_static("RANK", [](nb::handle) { return OpenMS::ProteinSchema::RANK; })
        .def_prop_ro_static("COVERAGE", [](nb::handle) { return OpenMS::ProteinSchema::COVERAGE; })
        .def_prop_ro_static("SEQUENCE", [](nb::handle) { return OpenMS::ProteinSchema::SEQUENCE; })
        .def_prop_ro_static("DESCRIPTION", [](nb::handle) { return OpenMS::ProteinSchema::DESCRIPTION; })
        .def_prop_ro_static("IS_DECOY", [](nb::handle) { return OpenMS::ProteinSchema::IS_DECOY; })
        .def_prop_ro_static("RUN_IDENTIFIER", [](nb::handle) { return OpenMS::ProteinSchema::RUN_IDENTIFIER; })
        .def_prop_ro_static("MODIFICATIONS", [](nb::handle) { return OpenMS::ProteinSchema::MODIFICATIONS; })
        .def_prop_ro_static("METAVALUES", [](nb::handle) { return OpenMS::ProteinSchema::METAVALUES; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::ProteinSchema::schema()); });

    nb::class_<OpenMS::ProteinGroupSchema>(m, "ProteinGroupSchema")
        .def_prop_ro_static("GROUP_TYPE", [](nb::handle) { return OpenMS::ProteinGroupSchema::GROUP_TYPE; })
        .def_prop_ro_static("PROBABILITY", [](nb::handle) { return OpenMS::ProteinGroupSchema::PROBABILITY; })
        .def_prop_ro_static("ACCESSIONS", [](nb::handle) { return OpenMS::ProteinGroupSchema::ACCESSIONS; })
        .def_prop_ro_static("RUN_IDENTIFIER", [](nb::handle) { return OpenMS::ProteinGroupSchema::RUN_IDENTIFIER; })
        .def_prop_ro_static("GROUP_INDEX", [](nb::handle) { return OpenMS::ProteinGroupSchema::GROUP_INDEX; })
        .def_prop_ro_static("FLOAT_DATA", [](nb::handle) { return OpenMS::ProteinGroupSchema::FLOAT_DATA; })
        .def_prop_ro_static("STRING_DATA", [](nb::handle) { return OpenMS::ProteinGroupSchema::STRING_DATA; })
        .def_prop_ro_static("INTEGER_DATA", [](nb::handle) { return OpenMS::ProteinGroupSchema::INTEGER_DATA; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::ProteinGroupSchema::schema()); });

    nb::class_<OpenMS::SearchParamsSchema>(m, "SearchParamsSchema")
        .def_prop_ro_static("RUN_IDENTIFIER", [](nb::handle) { return OpenMS::SearchParamsSchema::RUN_IDENTIFIER; })
        .def_prop_ro_static("SEARCH_ENGINE", [](nb::handle) { return OpenMS::SearchParamsSchema::SEARCH_ENGINE; })
        .def_prop_ro_static("SEARCH_ENGINE_VERSION", [](nb::handle) { return OpenMS::SearchParamsSchema::SEARCH_ENGINE_VERSION; })
        .def_prop_ro_static("INFERENCE_ENGINE", [](nb::handle) { return OpenMS::SearchParamsSchema::INFERENCE_ENGINE; })
        .def_prop_ro_static("INFERENCE_ENGINE_VERSION", [](nb::handle) { return OpenMS::SearchParamsSchema::INFERENCE_ENGINE_VERSION; })
        .def_prop_ro_static("DATE", [](nb::handle) { return OpenMS::SearchParamsSchema::DATE; })
        .def_prop_ro_static("SCORE_TYPE", [](nb::handle) { return OpenMS::SearchParamsSchema::SCORE_TYPE; })
        .def_prop_ro_static("HIGHER_SCORE_BETTER", [](nb::handle) { return OpenMS::SearchParamsSchema::HIGHER_SCORE_BETTER; })
        .def_prop_ro_static("SIGNIFICANCE_THRESHOLD", [](nb::handle) { return OpenMS::SearchParamsSchema::SIGNIFICANCE_THRESHOLD; })
        .def_prop_ro_static("DB", [](nb::handle) { return OpenMS::SearchParamsSchema::DB; })
        .def_prop_ro_static("DB_VERSION", [](nb::handle) { return OpenMS::SearchParamsSchema::DB_VERSION; })
        .def_prop_ro_static("TAXONOMY", [](nb::handle) { return OpenMS::SearchParamsSchema::TAXONOMY; })
        .def_prop_ro_static("CHARGES", [](nb::handle) { return OpenMS::SearchParamsSchema::CHARGES; })
        .def_prop_ro_static("MASS_TYPE", [](nb::handle) { return OpenMS::SearchParamsSchema::MASS_TYPE; })
        .def_prop_ro_static("PRECURSOR_MASS_TOLERANCE", [](nb::handle) { return OpenMS::SearchParamsSchema::PRECURSOR_MASS_TOLERANCE; })
        .def_prop_ro_static("PRECURSOR_MASS_TOLERANCE_PPM", [](nb::handle) { return OpenMS::SearchParamsSchema::PRECURSOR_MASS_TOLERANCE_PPM; })
        .def_prop_ro_static("FRAGMENT_MASS_TOLERANCE", [](nb::handle) { return OpenMS::SearchParamsSchema::FRAGMENT_MASS_TOLERANCE; })
        .def_prop_ro_static("FRAGMENT_MASS_TOLERANCE_PPM", [](nb::handle) { return OpenMS::SearchParamsSchema::FRAGMENT_MASS_TOLERANCE_PPM; })
        .def_prop_ro_static("DIGESTION_ENZYME", [](nb::handle) { return OpenMS::SearchParamsSchema::DIGESTION_ENZYME; })
        .def_prop_ro_static("ENZYME_TERM_SPECIFICITY", [](nb::handle) { return OpenMS::SearchParamsSchema::ENZYME_TERM_SPECIFICITY; })
        .def_prop_ro_static("MISSED_CLEAVAGES", [](nb::handle) { return OpenMS::SearchParamsSchema::MISSED_CLEAVAGES; })
        .def_prop_ro_static("FIXED_MODIFICATIONS", [](nb::handle) { return OpenMS::SearchParamsSchema::FIXED_MODIFICATIONS; })
        .def_prop_ro_static("VARIABLE_MODIFICATIONS", [](nb::handle) { return OpenMS::SearchParamsSchema::VARIABLE_MODIFICATIONS; })
        .def_prop_ro_static("PRIMARY_MS_RUN_PATHS", [](nb::handle) { return OpenMS::SearchParamsSchema::PRIMARY_MS_RUN_PATHS; })
        .def_prop_ro_static("METAVALUES", [](nb::handle) { return OpenMS::SearchParamsSchema::METAVALUES; })
        .def_prop_ro_static("SP_METAVALUES", [](nb::handle) { return OpenMS::SearchParamsSchema::SP_METAVALUES; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::SearchParamsSchema::schema()); });

    nb::class_<OpenMS::FeatureSchema>(m, "FeatureSchema")
        .def_prop_ro_static("UNIQUE_ID", [](nb::handle) { return OpenMS::FeatureSchema::UNIQUE_ID; })
        .def_prop_ro_static("PARENT_FEATURE_ID", [](nb::handle) { return OpenMS::FeatureSchema::PARENT_FEATURE_ID; })
        .def_prop_ro_static("DEPTH", [](nb::handle) { return OpenMS::FeatureSchema::DEPTH; })
        .def_prop_ro_static("RT", [](nb::handle) { return OpenMS::FeatureSchema::RT; })
        .def_prop_ro_static("MZ", [](nb::handle) { return OpenMS::FeatureSchema::MZ; })
        .def_prop_ro_static("INTENSITY", [](nb::handle) { return OpenMS::FeatureSchema::INTENSITY; })
        .def_prop_ro_static("CHARGE", [](nb::handle) { return OpenMS::FeatureSchema::CHARGE; })
        .def_prop_ro_static("QUALITY", [](nb::handle) { return OpenMS::FeatureSchema::QUALITY; })
        .def_prop_ro_static("QUALITY_RT", [](nb::handle) { return OpenMS::FeatureSchema::QUALITY_RT; })
        .def_prop_ro_static("QUALITY_MZ", [](nb::handle) { return OpenMS::FeatureSchema::QUALITY_MZ; })
        .def_prop_ro_static("WIDTH", [](nb::handle) { return OpenMS::FeatureSchema::WIDTH; })
        .def_prop_ro_static("RT_BB_MIN", [](nb::handle) { return OpenMS::FeatureSchema::RT_BB_MIN; })
        .def_prop_ro_static("RT_BB_MAX", [](nb::handle) { return OpenMS::FeatureSchema::RT_BB_MAX; })
        .def_prop_ro_static("MZ_BB_MIN", [](nb::handle) { return OpenMS::FeatureSchema::MZ_BB_MIN; })
        .def_prop_ro_static("MZ_BB_MAX", [](nb::handle) { return OpenMS::FeatureSchema::MZ_BB_MAX; })
        .def_prop_ro_static("CONVEX_HULLS", [](nb::handle) { return OpenMS::FeatureSchema::CONVEX_HULLS; })
        .def_prop_ro_static("METAVALUES", [](nb::handle) { return OpenMS::FeatureSchema::METAVALUES; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::FeatureSchema::schema()); });

    nb::class_<OpenMS::ConsensusFeatureSchema>(m, "ConsensusFeatureSchema")
        .def_prop_ro_static("UNIQUE_ID", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::UNIQUE_ID; })
        .def_prop_ro_static("RT", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::RT; })
        .def_prop_ro_static("MZ", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::MZ; })
        .def_prop_ro_static("INTENSITY", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::INTENSITY; })
        .def_prop_ro_static("CHARGE", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::CHARGE; })
        .def_prop_ro_static("QUALITY", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::QUALITY; })
        .def_prop_ro_static("WIDTH", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::WIDTH; })
        .def_prop_ro_static("HANDLES", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::HANDLES; })
        .def_prop_ro_static("METAVALUES", [](nb::handle) { return OpenMS::ConsensusFeatureSchema::METAVALUES; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::ConsensusFeatureSchema::schema()); });

    nb::class_<OpenMS::PSMSchema>(m, "PSMSchema")
        .def_prop_ro_static("SEQUENCE", [](nb::handle) { return OpenMS::PSMSchema::SEQUENCE; })
        .def_prop_ro_static("PEPTIDOFORM", [](nb::handle) { return OpenMS::PSMSchema::PEPTIDOFORM; })
        .def_prop_ro_static("MODIFICATIONS", [](nb::handle) { return OpenMS::PSMSchema::MODIFICATIONS; })
        .def_prop_ro_static("PRECURSOR_CHARGE", [](nb::handle) { return OpenMS::PSMSchema::PRECURSOR_CHARGE; })
        .def_prop_ro_static("POSTERIOR_ERROR_PROBABILITY", [](nb::handle) { return OpenMS::PSMSchema::POSTERIOR_ERROR_PROBABILITY; })
        .def_prop_ro_static("IS_DECOY", [](nb::handle) { return OpenMS::PSMSchema::IS_DECOY; })
        .def_prop_ro_static("CALCULATED_MZ", [](nb::handle) { return OpenMS::PSMSchema::CALCULATED_MZ; })
        .def_prop_ro_static("OBSERVED_MZ", [](nb::handle) { return OpenMS::PSMSchema::OBSERVED_MZ; })
        .def_prop_ro_static("ADDITIONAL_SCORES", [](nb::handle) { return OpenMS::PSMSchema::ADDITIONAL_SCORES; })
        .def_prop_ro_static("PROTEIN_ACCESSIONS", [](nb::handle) { return OpenMS::PSMSchema::PROTEIN_ACCESSIONS; })
        .def_prop_ro_static("PREDICTED_RT", [](nb::handle) { return OpenMS::PSMSchema::PREDICTED_RT; })
        .def_prop_ro_static("REFERENCE_FILE_NAME", [](nb::handle) { return OpenMS::PSMSchema::REFERENCE_FILE_NAME; })
        .def_prop_ro_static("CV_PARAMS", [](nb::handle) { return OpenMS::PSMSchema::CV_PARAMS; })
        .def_prop_ro_static("SCAN", [](nb::handle) { return OpenMS::PSMSchema::SCAN; })
        .def_prop_ro_static("RT", [](nb::handle) { return OpenMS::PSMSchema::RT; })
        .def_prop_ro_static("ION_MOBILITY", [](nb::handle) { return OpenMS::PSMSchema::ION_MOBILITY; })
        .def_prop_ro_static("SPECTRUM_REFERENCE", [](nb::handle) { return OpenMS::PSMSchema::SPECTRUM_REFERENCE; })
        .def_prop_ro_static("SCORE", [](nb::handle) { return OpenMS::PSMSchema::SCORE; })
        .def_prop_ro_static("SCORE_TYPE", [](nb::handle) { return OpenMS::PSMSchema::SCORE_TYPE; })
        .def_prop_ro_static("HIGHER_SCORE_BETTER", [](nb::handle) { return OpenMS::PSMSchema::HIGHER_SCORE_BETTER; })
        .def_prop_ro_static("RANK", [](nb::handle) { return OpenMS::PSMSchema::RANK; })
        .def_prop_ro_static("PEPTIDE_IDENTIFICATION_INDEX", [](nb::handle) { return OpenMS::PSMSchema::PEPTIDE_IDENTIFICATION_INDEX; })
        .def_prop_ro_static("PSM_METAVALUES", [](nb::handle) { return OpenMS::PSMSchema::PSM_METAVALUES; })
        .def_prop_ro_static("SPECTRUM_METAVALUES", [](nb::handle) { return OpenMS::PSMSchema::SPECTRUM_METAVALUES; })
        .def_prop_ro_static("RUN_IDENTIFIER", [](nb::handle) { return OpenMS::PSMSchema::RUN_IDENTIFIER; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::PSMSchema::schema()); });

    nb::class_<OpenMS::ConsensusFeatureExportSchema>(m, "ConsensusFeatureExportSchema")
        .def_prop_ro_static("SEQUENCE", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::SEQUENCE; })
        .def_prop_ro_static("PEPTIDOFORM", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::PEPTIDOFORM; })
        .def_prop_ro_static("MODIFICATIONS", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::MODIFICATIONS; })
        .def_prop_ro_static("PRECURSOR_CHARGE", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::PRECURSOR_CHARGE; })
        .def_prop_ro_static("CALCULATED_MZ", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::CALCULATED_MZ; })
        .def_prop_ro_static("OBSERVED_MZ", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::OBSERVED_MZ; })
        .def_prop_ro_static("RT", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::RT; })
        .def_prop_ro_static("POSTERIOR_ERROR_PROBABILITY", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::POSTERIOR_ERROR_PROBABILITY; })
        .def_prop_ro_static("IS_DECOY", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::IS_DECOY; })
        .def_prop_ro_static("ADDITIONAL_SCORES", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::ADDITIONAL_SCORES; })
        .def_prop_ro_static("PREDICTED_RT", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::PREDICTED_RT; })
        .def_prop_ro_static("REFERENCE_FILE_NAME", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::REFERENCE_FILE_NAME; })
        .def_prop_ro_static("CV_PARAMS", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::CV_PARAMS; })
        .def_prop_ro_static("SCAN", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::SCAN; })
        .def_prop_ro_static("ION_MOBILITY", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::ION_MOBILITY; })
        .def_prop_ro_static("START_ION_MOBILITY", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::START_ION_MOBILITY; })
        .def_prop_ro_static("STOP_ION_MOBILITY", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::STOP_ION_MOBILITY; })
        .def_prop_ro_static("INTENSITIES", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::INTENSITIES; })
        .def_prop_ro_static("ADDITIONAL_INTENSITIES", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::ADDITIONAL_INTENSITIES; })
        .def_prop_ro_static("PG_ACCESSIONS", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::PG_ACCESSIONS; })
        .def_prop_ro_static("ANCHOR_PROTEIN", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::ANCHOR_PROTEIN; })
        .def_prop_ro_static("UNIQUE", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::UNIQUE; })
        .def_prop_ro_static("PG_GLOBAL_QVALUE", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::PG_GLOBAL_QVALUE; })
        .def_prop_ro_static("GG_ACCESSIONS", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::GG_ACCESSIONS; })
        .def_prop_ro_static("GG_NAMES", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::GG_NAMES; })
        .def_prop_ro_static("SCAN_REFERENCE_FILE_NAME", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::SCAN_REFERENCE_FILE_NAME; })
        .def_prop_ro_static("RT_START", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::RT_START; })
        .def_prop_ro_static("RT_STOP", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::RT_STOP; })
        .def_prop_ro_static("QUALITY", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::QUALITY; })
        .def_prop_ro_static("SCORE", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::SCORE; })
        .def_prop_ro_static("SCORE_TYPE", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::SCORE_TYPE; })
        .def_prop_ro_static("SPECTRUM_REFERENCE", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::SPECTRUM_REFERENCE; })
        .def_prop_ro_static("FEATURE_METAVALUES", [](nb::handle) { return OpenMS::ConsensusFeatureExportSchema::FEATURE_METAVALUES; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::ConsensusFeatureExportSchema::schema()); });

    nb::class_<OpenMS::SpectraLongSchema>(m, "SpectraLongSchema")
        .def_prop_ro_static("MZ", [](nb::handle) { return OpenMS::SpectraLongSchema::MZ; })
        .def_prop_ro_static("INTENSITY", [](nb::handle) { return OpenMS::SpectraLongSchema::INTENSITY; })
        .def_prop_ro_static("RT", [](nb::handle) { return OpenMS::SpectraLongSchema::RT; })
        .def_prop_ro_static("ION_MOBILITY", [](nb::handle) { return OpenMS::SpectraLongSchema::ION_MOBILITY; })
        .def_prop_ro_static("SPECTRUM_INDEX", [](nb::handle) { return OpenMS::SpectraLongSchema::SPECTRUM_INDEX; })
        .def_prop_ro_static("MS_LEVEL", [](nb::handle) { return OpenMS::SpectraLongSchema::MS_LEVEL; })
        .def_prop_ro_static("NATIVE_ID", [](nb::handle) { return OpenMS::SpectraLongSchema::NATIVE_ID; })
        .def_prop_ro_static("PRECURSOR_MZ", [](nb::handle) { return OpenMS::SpectraLongSchema::PRECURSOR_MZ; })
        .def_prop_ro_static("PRECURSOR_CHARGE", [](nb::handle) { return OpenMS::SpectraLongSchema::PRECURSOR_CHARGE; })
        .def_prop_ro_static("PRECURSOR_INTENSITY", [](nb::handle) { return OpenMS::SpectraLongSchema::PRECURSOR_INTENSITY; })
        .def_prop_ro_static("ISOLATION_LOWER", [](nb::handle) { return OpenMS::SpectraLongSchema::ISOLATION_LOWER; })
        .def_prop_ro_static("ISOLATION_UPPER", [](nb::handle) { return OpenMS::SpectraLongSchema::ISOLATION_UPPER; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::SpectraLongSchema::schema()); });

    nb::class_<OpenMS::SpectraSemiWideSchema>(m, "SpectraSemiWideSchema")
        .def_prop_ro_static("SPECTRUM_INDEX", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::SPECTRUM_INDEX; })
        .def_prop_ro_static("RT", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::RT; })
        .def_prop_ro_static("MS_LEVEL", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::MS_LEVEL; })
        .def_prop_ro_static("NATIVE_ID", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::NATIVE_ID; })
        .def_prop_ro_static("MZ", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::MZ; })
        .def_prop_ro_static("INTENSITY", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::INTENSITY; })
        .def_prop_ro_static("ION_MOBILITY", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::ION_MOBILITY; })
        .def_prop_ro_static("PRECURSOR_MZ", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::PRECURSOR_MZ; })
        .def_prop_ro_static("PRECURSOR_CHARGE", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::PRECURSOR_CHARGE; })
        .def_prop_ro_static("PRECURSOR_INTENSITY", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::PRECURSOR_INTENSITY; })
        .def_prop_ro_static("ISOLATION_LOWER", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::ISOLATION_LOWER; })
        .def_prop_ro_static("ISOLATION_UPPER", [](nb::handle) { return OpenMS::SpectraSemiWideSchema::ISOLATION_UPPER; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::SpectraSemiWideSchema::schema()); });

    nb::class_<OpenMS::ChromatogramSchema>(m, "ChromatogramSchema")
        .def_prop_ro_static("RT", [](nb::handle) { return OpenMS::ChromatogramSchema::RT; })
        .def_prop_ro_static("INTENSITY", [](nb::handle) { return OpenMS::ChromatogramSchema::INTENSITY; })
        .def_prop_ro_static("CHROMATOGRAM_INDEX", [](nb::handle) { return OpenMS::ChromatogramSchema::CHROMATOGRAM_INDEX; })
        .def_prop_ro_static("NATIVE_ID", [](nb::handle) { return OpenMS::ChromatogramSchema::NATIVE_ID; })
        .def_prop_ro_static("PRECURSOR_MZ", [](nb::handle) { return OpenMS::ChromatogramSchema::PRECURSOR_MZ; })
        .def_prop_ro_static("PRODUCT_MZ", [](nb::handle) { return OpenMS::ChromatogramSchema::PRODUCT_MZ; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::ChromatogramSchema::schema()); });

    nb::class_<OpenMS::OSWPrecursorSchema>(m, "OSWPrecursorSchema")
        .def_prop_ro_static("PRECURSOR_ID", [](nb::handle) { return OpenMS::OSWPrecursorSchema::PRECURSOR_ID; })
        .def_prop_ro_static("PRECURSOR_MZ", [](nb::handle) { return OpenMS::OSWPrecursorSchema::PRECURSOR_MZ; })
        .def_prop_ro_static("CHARGE", [](nb::handle) { return OpenMS::OSWPrecursorSchema::CHARGE; })
        .def_prop_ro_static("LIBRARY_RT", [](nb::handle) { return OpenMS::OSWPrecursorSchema::LIBRARY_RT; })
        .def_prop_ro_static("LIBRARY_DRIFT_TIME", [](nb::handle) { return OpenMS::OSWPrecursorSchema::LIBRARY_DRIFT_TIME; })
        .def_prop_ro_static("DECOY", [](nb::handle) { return OpenMS::OSWPrecursorSchema::DECOY; })
        .def_prop_ro_static("TRAML_ID", [](nb::handle) { return OpenMS::OSWPrecursorSchema::TRAML_ID; })
        .def_prop_ro_static("MODIFIED_SEQUENCE", [](nb::handle) { return OpenMS::OSWPrecursorSchema::MODIFIED_SEQUENCE; })
        .def_prop_ro_static("UNMODIFIED_SEQUENCE", [](nb::handle) { return OpenMS::OSWPrecursorSchema::UNMODIFIED_SEQUENCE; })
        .def_prop_ro_static("PROTEIN_ACCESSIONS", [](nb::handle) { return OpenMS::OSWPrecursorSchema::PROTEIN_ACCESSIONS; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::OSWPrecursorSchema::schema()); });

    nb::class_<OpenMS::OSWTransitionSchema>(m, "OSWTransitionSchema")
        .def_prop_ro_static("TRANSITION_ID", [](nb::handle) { return OpenMS::OSWTransitionSchema::TRANSITION_ID; })
        .def_prop_ro_static("PRECURSOR_ID", [](nb::handle) { return OpenMS::OSWTransitionSchema::PRECURSOR_ID; })
        .def_prop_ro_static("TRAML_ID", [](nb::handle) { return OpenMS::OSWTransitionSchema::TRAML_ID; })
        .def_prop_ro_static("PRODUCT_MZ", [](nb::handle) { return OpenMS::OSWTransitionSchema::PRODUCT_MZ; })
        .def_prop_ro_static("CHARGE", [](nb::handle) { return OpenMS::OSWTransitionSchema::CHARGE; })
        .def_prop_ro_static("TYPE", [](nb::handle) { return OpenMS::OSWTransitionSchema::TYPE; })
        .def_prop_ro_static("ANNOTATION", [](nb::handle) { return OpenMS::OSWTransitionSchema::ANNOTATION; })
        .def_prop_ro_static("ORDINAL", [](nb::handle) { return OpenMS::OSWTransitionSchema::ORDINAL; })
        .def_prop_ro_static("DETECTING", [](nb::handle) { return OpenMS::OSWTransitionSchema::DETECTING; })
        .def_prop_ro_static("IDENTIFYING", [](nb::handle) { return OpenMS::OSWTransitionSchema::IDENTIFYING; })
        .def_prop_ro_static("QUANTIFYING", [](nb::handle) { return OpenMS::OSWTransitionSchema::QUANTIFYING; })
        .def_prop_ro_static("LIBRARY_INTENSITY", [](nb::handle) { return OpenMS::OSWTransitionSchema::LIBRARY_INTENSITY; })
        .def_prop_ro_static("DECOY", [](nb::handle) { return OpenMS::OSWTransitionSchema::DECOY; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::OSWTransitionSchema::schema()); });

    nb::class_<OpenMS::OSWFeaturePrecursorSchema>(m, "OSWFeaturePrecursorSchema")
        .def_prop_ro_static("FEATURE_ID", [](nb::handle) { return OpenMS::OSWFeaturePrecursorSchema::FEATURE_ID; })
        .def_prop_ro_static("RUN_ID", [](nb::handle) { return OpenMS::OSWFeaturePrecursorSchema::RUN_ID; })
        .def_prop_ro_static("PRECURSOR_ISOTOPE", [](nb::handle) { return OpenMS::OSWFeaturePrecursorSchema::PRECURSOR_ISOTOPE; })
        .def_prop_ro_static("PRECURSOR_AREA_INTENSITY", [](nb::handle) { return OpenMS::OSWFeaturePrecursorSchema::PRECURSOR_AREA_INTENSITY; })
        .def_prop_ro_static("PRECURSOR_APEX_INTENSITY", [](nb::handle) { return OpenMS::OSWFeaturePrecursorSchema::PRECURSOR_APEX_INTENSITY; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::OSWFeaturePrecursorSchema::schema()); });

    nb::class_<OpenMS::OSWRunSchema>(m, "OSWRunSchema")
        .def_prop_ro_static("RUN_ID", [](nb::handle) { return OpenMS::OSWRunSchema::RUN_ID; })
        .def_prop_ro_static("FILENAME", [](nb::handle) { return OpenMS::OSWRunSchema::FILENAME; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::OSWRunSchema::schema()); });

    nb::class_<OpenMS::OSWFeatureSchema>(m, "OSWFeatureSchema")
        .def_prop_ro_static("FEATURE_ID", [](nb::handle) { return OpenMS::OSWFeatureSchema::FEATURE_ID; })
        .def_prop_ro_static("RUN_ID", [](nb::handle) { return OpenMS::OSWFeatureSchema::RUN_ID; })
        .def_prop_ro_static("PRECURSOR_ID", [](nb::handle) { return OpenMS::OSWFeatureSchema::PRECURSOR_ID; })
        .def_prop_ro_static("EXP_RT", [](nb::handle) { return OpenMS::OSWFeatureSchema::EXP_RT; })
        .def_prop_ro_static("EXP_IM", [](nb::handle) { return OpenMS::OSWFeatureSchema::EXP_IM; })
        .def_prop_ro_static("NORM_RT", [](nb::handle) { return OpenMS::OSWFeatureSchema::NORM_RT; })
        .def_prop_ro_static("DELTA_RT", [](nb::handle) { return OpenMS::OSWFeatureSchema::DELTA_RT; })
        .def_prop_ro_static("LEFT_WIDTH", [](nb::handle) { return OpenMS::OSWFeatureSchema::LEFT_WIDTH; })
        .def_prop_ro_static("RIGHT_WIDTH", [](nb::handle) { return OpenMS::OSWFeatureSchema::RIGHT_WIDTH; })
        .def_prop_ro_static("EXP_IM_LEFTWIDTH", [](nb::handle) { return OpenMS::OSWFeatureSchema::EXP_IM_LEFTWIDTH; })
        .def_prop_ro_static("EXP_IM_RIGHTWIDTH", [](nb::handle) { return OpenMS::OSWFeatureSchema::EXP_IM_RIGHTWIDTH; })
        .def_prop_ro_static("MS1_AREA_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS1_AREA_INTENSITY; })
        .def_prop_ro_static("MS1_APEX_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS1_APEX_INTENSITY; })
        .def_prop_ro_static("MS1_EXP_IM", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS1_EXP_IM; })
        .def_prop_ro_static("MS1_DELTA_IM", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS1_DELTA_IM; })
        .def_prop_ro_static("VAR_MS1_MASSDEV_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_MASSDEV_SCORE; })
        .def_prop_ro_static("VAR_MS1_IM_MS1_DELTA_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_IM_MS1_DELTA_SCORE; })
        .def_prop_ro_static("VAR_MS1_MI_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_MI_SCORE; })
        .def_prop_ro_static("VAR_MS1_MI_CONTRAST_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_MI_CONTRAST_SCORE; })
        .def_prop_ro_static("VAR_MS1_MI_COMBINED_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_MI_COMBINED_SCORE; })
        .def_prop_ro_static("VAR_MS1_ISOTOPE_CORRELATION_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_ISOTOPE_CORRELATION_SCORE; })
        .def_prop_ro_static("VAR_MS1_ISOTOPE_OVERLAP_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_ISOTOPE_OVERLAP_SCORE; })
        .def_prop_ro_static("VAR_MS1_XCORR_COELUTION", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_XCORR_COELUTION; })
        .def_prop_ro_static("VAR_MS1_XCORR_COELUTION_CONTRAST", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_XCORR_COELUTION_CONTRAST; })
        .def_prop_ro_static("VAR_MS1_XCORR_COELUTION_COMBINED", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_XCORR_COELUTION_COMBINED; })
        .def_prop_ro_static("VAR_MS1_XCORR_SHAPE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_XCORR_SHAPE; })
        .def_prop_ro_static("VAR_MS1_XCORR_SHAPE_CONTRAST", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_XCORR_SHAPE_CONTRAST; })
        .def_prop_ro_static("VAR_MS1_XCORR_SHAPE_COMBINED", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS1_XCORR_SHAPE_COMBINED; })
        .def_prop_ro_static("MS2_AREA_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS2_AREA_INTENSITY; })
        .def_prop_ro_static("MS2_TOTAL_AREA_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS2_TOTAL_AREA_INTENSITY; })
        .def_prop_ro_static("MS2_APEX_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS2_APEX_INTENSITY; })
        .def_prop_ro_static("MS2_EXP_IM", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS2_EXP_IM; })
        .def_prop_ro_static("MS2_EXP_IM_LEFTWIDTH", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS2_EXP_IM_LEFTWIDTH; })
        .def_prop_ro_static("MS2_EXP_IM_RIGHTWIDTH", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS2_EXP_IM_RIGHTWIDTH; })
        .def_prop_ro_static("MS2_DELTA_IM", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS2_DELTA_IM; })
        .def_prop_ro_static("MS2_TOTAL_MI", [](nb::handle) { return OpenMS::OSWFeatureSchema::MS2_TOTAL_MI; })
        .def_prop_ro_static("VAR_MS2_BSERIES_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_BSERIES_SCORE; })
        .def_prop_ro_static("VAR_MS2_DOTPROD_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_DOTPROD_SCORE; })
        .def_prop_ro_static("VAR_MS2_INTENSITY_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_INTENSITY_SCORE; })
        .def_prop_ro_static("VAR_MS2_ISOTOPE_CORRELATION_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_ISOTOPE_CORRELATION_SCORE; })
        .def_prop_ro_static("VAR_MS2_ISOTOPE_OVERLAP_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_ISOTOPE_OVERLAP_SCORE; })
        .def_prop_ro_static("VAR_MS2_LIBRARY_CORR", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_LIBRARY_CORR; })
        .def_prop_ro_static("VAR_MS2_LIBRARY_DOTPROD", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_LIBRARY_DOTPROD; })
        .def_prop_ro_static("VAR_MS2_LIBRARY_MANHATTAN", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_LIBRARY_MANHATTAN; })
        .def_prop_ro_static("VAR_MS2_LIBRARY_RMSD", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_LIBRARY_RMSD; })
        .def_prop_ro_static("VAR_MS2_LIBRARY_ROOTMEANSQUARE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_LIBRARY_ROOTMEANSQUARE; })
        .def_prop_ro_static("VAR_MS2_LIBRARY_SANGLE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_LIBRARY_SANGLE; })
        .def_prop_ro_static("VAR_MS2_LOG_SN_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_LOG_SN_SCORE; })
        .def_prop_ro_static("VAR_MS2_MANHATTAN_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_MANHATTAN_SCORE; })
        .def_prop_ro_static("VAR_MS2_MASSDEV_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_MASSDEV_SCORE; })
        .def_prop_ro_static("VAR_MS2_MASSDEV_SCORE_WEIGHTED", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_MASSDEV_SCORE_WEIGHTED; })
        .def_prop_ro_static("VAR_MS2_MI_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_MI_SCORE; })
        .def_prop_ro_static("VAR_MS2_MI_WEIGHTED_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_MI_WEIGHTED_SCORE; })
        .def_prop_ro_static("VAR_MS2_MI_RATIO_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_MI_RATIO_SCORE; })
        .def_prop_ro_static("VAR_MS2_NORM_RT_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_NORM_RT_SCORE; })
        .def_prop_ro_static("VAR_MS2_XCORR_COELUTION", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_XCORR_COELUTION; })
        .def_prop_ro_static("VAR_MS2_XCORR_COELUTION_WEIGHTED", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_XCORR_COELUTION_WEIGHTED; })
        .def_prop_ro_static("VAR_MS2_XCORR_SHAPE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_XCORR_SHAPE; })
        .def_prop_ro_static("VAR_MS2_XCORR_SHAPE_WEIGHTED", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_XCORR_SHAPE_WEIGHTED; })
        .def_prop_ro_static("VAR_MS2_YSERIES_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_YSERIES_SCORE; })
        .def_prop_ro_static("VAR_MS2_ELUTION_MODEL_FIT_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_ELUTION_MODEL_FIT_SCORE; })
        .def_prop_ro_static("VAR_MS2_IM_XCORR_SHAPE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_IM_XCORR_SHAPE; })
        .def_prop_ro_static("VAR_MS2_IM_XCORR_COELUTION", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_IM_XCORR_COELUTION; })
        .def_prop_ro_static("VAR_MS2_IM_DELTA_SCORE", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_IM_DELTA_SCORE; })
        .def_prop_ro_static("VAR_MS2_IM_LOG_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureSchema::VAR_MS2_IM_LOG_INTENSITY; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::OSWFeatureSchema::schema()); });

    nb::class_<OpenMS::OSWFeatureTransitionSchema>(m, "OSWFeatureTransitionSchema")
        .def_prop_ro_static("FEATURE_ID", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::FEATURE_ID; })
        .def_prop_ro_static("RUN_ID", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::RUN_ID; })
        .def_prop_ro_static("TRANSITION_ID", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::TRANSITION_ID; })
        .def_prop_ro_static("AREA_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::AREA_INTENSITY; })
        .def_prop_ro_static("TOTAL_AREA_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::TOTAL_AREA_INTENSITY; })
        .def_prop_ro_static("APEX_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::APEX_INTENSITY; })
        .def_prop_ro_static("APEX_RT", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::APEX_RT; })
        .def_prop_ro_static("RT_FWHM", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::RT_FWHM; })
        .def_prop_ro_static("MASSERROR_PPM", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::MASSERROR_PPM; })
        .def_prop_ro_static("TOTAL_MI", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::TOTAL_MI; })
        .def_prop_ro_static("VAR_INTENSITY_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_INTENSITY_SCORE; })
        .def_prop_ro_static("VAR_INTENSITY_RATIO_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_INTENSITY_RATIO_SCORE; })
        .def_prop_ro_static("VAR_LOG_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_LOG_INTENSITY; })
        .def_prop_ro_static("VAR_XCORR_COELUTION", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_XCORR_COELUTION; })
        .def_prop_ro_static("VAR_XCORR_SHAPE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_XCORR_SHAPE; })
        .def_prop_ro_static("VAR_LOG_SN_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_LOG_SN_SCORE; })
        .def_prop_ro_static("VAR_MASSDEV_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_MASSDEV_SCORE; })
        .def_prop_ro_static("VAR_MI_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_MI_SCORE; })
        .def_prop_ro_static("VAR_MI_RATIO_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_MI_RATIO_SCORE; })
        .def_prop_ro_static("VAR_ISOTOPE_CORRELATION_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_ISOTOPE_CORRELATION_SCORE; })
        .def_prop_ro_static("VAR_ISOTOPE_OVERLAP_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_ISOTOPE_OVERLAP_SCORE; })
        .def_prop_ro_static("EXP_IM", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::EXP_IM; })
        .def_prop_ro_static("EXP_IM_LEFTWIDTH", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::EXP_IM_LEFTWIDTH; })
        .def_prop_ro_static("EXP_IM_RIGHTWIDTH", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::EXP_IM_RIGHTWIDTH; })
        .def_prop_ro_static("DELTA_IM", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::DELTA_IM; })
        .def_prop_ro_static("VAR_IM_DELTA_SCORE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_IM_DELTA_SCORE; })
        .def_prop_ro_static("VAR_IM_LOG_INTENSITY", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_IM_LOG_INTENSITY; })
        .def_prop_ro_static("VAR_IM_XCORR_COELUTION_CONTRAST", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_IM_XCORR_COELUTION_CONTRAST; })
        .def_prop_ro_static("VAR_IM_XCORR_SHAPE_CONTRAST", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_IM_XCORR_SHAPE_CONTRAST; })
        .def_prop_ro_static("VAR_IM_XCORR_COELUTION_COMBINED", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_IM_XCORR_COELUTION_COMBINED; })
        .def_prop_ro_static("VAR_IM_XCORR_SHAPE_COMBINED", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::VAR_IM_XCORR_SHAPE_COMBINED; })
        .def_prop_ro_static("START_POSITION_AT_5", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::START_POSITION_AT_5; })
        .def_prop_ro_static("END_POSITION_AT_5", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::END_POSITION_AT_5; })
        .def_prop_ro_static("START_POSITION_AT_10", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::START_POSITION_AT_10; })
        .def_prop_ro_static("END_POSITION_AT_10", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::END_POSITION_AT_10; })
        .def_prop_ro_static("START_POSITION_AT_50", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::START_POSITION_AT_50; })
        .def_prop_ro_static("END_POSITION_AT_50", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::END_POSITION_AT_50; })
        .def_prop_ro_static("TOTAL_WIDTH", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::TOTAL_WIDTH; })
        .def_prop_ro_static("TAILING_FACTOR", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::TAILING_FACTOR; })
        .def_prop_ro_static("ASYMMETRY_FACTOR", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::ASYMMETRY_FACTOR; })
        .def_prop_ro_static("SLOPE_OF_BASELINE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::SLOPE_OF_BASELINE; })
        .def_prop_ro_static("BASELINE_DELTA_2_HEIGHT", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::BASELINE_DELTA_2_HEIGHT; })
        .def_prop_ro_static("POINTS_ACROSS_BASELINE", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::POINTS_ACROSS_BASELINE; })
        .def_prop_ro_static("POINTS_ACROSS_HALF_HEIGHT", [](nb::handle) { return OpenMS::OSWFeatureTransitionSchema::POINTS_ACROSS_HALF_HEIGHT; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::OSWFeatureTransitionSchema::schema()); });

    nb::class_<OpenMS::XICSchema>(m, "XICSchema")
        .def_prop_ro_static("RUN_ID", [](nb::handle) { return OpenMS::XICSchema::RUN_ID; })
        .def_prop_ro_static("SOURCE_FILE", [](nb::handle) { return OpenMS::XICSchema::SOURCE_FILE; })
        .def_prop_ro_static("MS_LEVEL", [](nb::handle) { return OpenMS::XICSchema::MS_LEVEL; })
        .def_prop_ro_static("PRECURSOR_ID", [](nb::handle) { return OpenMS::XICSchema::PRECURSOR_ID; })
        .def_prop_ro_static("TRANSITION_ID", [](nb::handle) { return OpenMS::XICSchema::TRANSITION_ID; })
        .def_prop_ro_static("MODIFIED_SEQUENCE", [](nb::handle) { return OpenMS::XICSchema::MODIFIED_SEQUENCE; })
        .def_prop_ro_static("PRECURSOR_CHARGE", [](nb::handle) { return OpenMS::XICSchema::PRECURSOR_CHARGE; })
        .def_prop_ro_static("PRODUCT_CHARGE", [](nb::handle) { return OpenMS::XICSchema::PRODUCT_CHARGE; })
        .def_prop_ro_static("DETECTING_TRANSITION", [](nb::handle) { return OpenMS::XICSchema::DETECTING_TRANSITION; })
        .def_prop_ro_static("PRECURSOR_DECOY", [](nb::handle) { return OpenMS::XICSchema::PRECURSOR_DECOY; })
        .def_prop_ro_static("PRODUCT_DECOY", [](nb::handle) { return OpenMS::XICSchema::PRODUCT_DECOY; })
        .def_prop_ro_static("TRANSITION_ORDINAL", [](nb::handle) { return OpenMS::XICSchema::TRANSITION_ORDINAL; })
        .def_prop_ro_static("TRANSITION_TYPE", [](nb::handle) { return OpenMS::XICSchema::TRANSITION_TYPE; })
        .def_prop_ro_static("ANNOTATION", [](nb::handle) { return OpenMS::XICSchema::ANNOTATION; })
        .def_prop_ro_static("RT_DATA", [](nb::handle) { return OpenMS::XICSchema::RT_DATA; })
        .def_prop_ro_static("INTENSITY_DATA", [](nb::handle) { return OpenMS::XICSchema::INTENSITY_DATA; })
        .def_prop_ro_static("RT_COMPRESSION", [](nb::handle) { return OpenMS::XICSchema::RT_COMPRESSION; })
        .def_prop_ro_static("INTENSITY_COMPRESSION", [](nb::handle) { return OpenMS::XICSchema::INTENSITY_COMPRESSION; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::XICSchema::schema()); });

    nb::class_<OpenMS::XIMSchema>(m, "XIMSchema")
        .def_prop_ro_static("RUN_ID", [](nb::handle) { return OpenMS::XIMSchema::RUN_ID; })
        .def_prop_ro_static("SOURCE_FILE", [](nb::handle) { return OpenMS::XIMSchema::SOURCE_FILE; })
        .def_prop_ro_static("MS_LEVEL", [](nb::handle) { return OpenMS::XIMSchema::MS_LEVEL; })
        .def_prop_ro_static("MOBILOGRAM_TYPE", [](nb::handle) { return OpenMS::XIMSchema::MOBILOGRAM_TYPE; })
        .def_prop_ro_static("PRECURSOR_ID", [](nb::handle) { return OpenMS::XIMSchema::PRECURSOR_ID; })
        .def_prop_ro_static("TRANSITION_ID", [](nb::handle) { return OpenMS::XIMSchema::TRANSITION_ID; })
        .def_prop_ro_static("FEATURE_ID", [](nb::handle) { return OpenMS::XIMSchema::FEATURE_ID; })
        .def_prop_ro_static("FEATURE_RT", [](nb::handle) { return OpenMS::XIMSchema::FEATURE_RT; })
        .def_prop_ro_static("MODIFIED_SEQUENCE", [](nb::handle) { return OpenMS::XIMSchema::MODIFIED_SEQUENCE; })
        .def_prop_ro_static("PRECURSOR_CHARGE", [](nb::handle) { return OpenMS::XIMSchema::PRECURSOR_CHARGE; })
        .def_prop_ro_static("PRODUCT_CHARGE", [](nb::handle) { return OpenMS::XIMSchema::PRODUCT_CHARGE; })
        .def_prop_ro_static("DETECTING_TRANSITION", [](nb::handle) { return OpenMS::XIMSchema::DETECTING_TRANSITION; })
        .def_prop_ro_static("PRECURSOR_DECOY", [](nb::handle) { return OpenMS::XIMSchema::PRECURSOR_DECOY; })
        .def_prop_ro_static("PRODUCT_DECOY", [](nb::handle) { return OpenMS::XIMSchema::PRODUCT_DECOY; })
        .def_prop_ro_static("TRANSITION_ORDINAL", [](nb::handle) { return OpenMS::XIMSchema::TRANSITION_ORDINAL; })
        .def_prop_ro_static("TRANSITION_TYPE", [](nb::handle) { return OpenMS::XIMSchema::TRANSITION_TYPE; })
        .def_prop_ro_static("ANNOTATION", [](nb::handle) { return OpenMS::XIMSchema::ANNOTATION; })
        .def_prop_ro_static("MOBILITY_DATA", [](nb::handle) { return OpenMS::XIMSchema::MOBILITY_DATA; })
        .def_prop_ro_static("INTENSITY_DATA", [](nb::handle) { return OpenMS::XIMSchema::INTENSITY_DATA; })
        .def_prop_ro_static("MOBILITY_COMPRESSION", [](nb::handle) { return OpenMS::XIMSchema::MOBILITY_COMPRESSION; })
        .def_prop_ro_static("INTENSITY_COMPRESSION", [](nb::handle) { return OpenMS::XIMSchema::INTENSITY_COMPRESSION; })
        .def_static("schema", []() { return schema_to_pyarrow(OpenMS::XIMSchema::schema()); });

    // -----------------------------------------------------------------------
    // ArrowSchemaValidation — ValidationResult and validate_arrow_schema
    // -----------------------------------------------------------------------

    nb::class_<OpenMS::ArrowSchemaValidation::ValidationResult>(m, "ValidationResult")
        .def_ro("valid", &OpenMS::ArrowSchemaValidation::ValidationResult::valid)
        .def_ro("errors", &OpenMS::ArrowSchemaValidation::ValidationResult::errors)
        .def("__str__", &OpenMS::ArrowSchemaValidation::ValidationResult::toString);

    m.def("validate_arrow_schema",
        [](nb::object py_table, nb::object py_schema, const std::string& mode) {
            auto table = pyarrow_to_table(py_table);
            auto schema = pyarrow_schema_to_arrow(py_schema);
            auto validation_mode = (mode == "subset")
                ? OpenMS::ArrowSchemaValidation::Mode::Subset
                : OpenMS::ArrowSchemaValidation::Mode::Strict;
            return OpenMS::ArrowSchemaValidation::validate(table, schema, validation_mode);
        },
        nb::arg("table"), nb::arg("schema"), nb::arg("mode") = "strict",
        "Validate an Arrow table against an expected schema.\n\n"
        "Parameters\n"
        "----------\n"
        "table : pyarrow.Table\n"
        "    The table to validate.\n"
        "schema : pyarrow.Schema\n"
        "    The expected schema.\n"
        "mode : str, optional\n"
        "    Validation mode: 'strict' (default) requires exact match,\n"
        "    'subset' only requires that table columns are a subset of the schema.\n\n"
        "Returns\n"
        "-------\n"
        "ValidationResult\n"
        "    Result with .valid (bool) and .errors (list of str).\n\n"
        ".. warning::\n"
        "    **EXPERIMENTAL API**: This function is experimental and may change."
    );
}
