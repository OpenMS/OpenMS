# Part of pyOpenMS. Copyright 2002-present, OpenMS Inc.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from Types cimport *
from ProteinIdentification cimport *
from MSExperimentArrowExport cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

# OpenMS ProteinIdentificationArrowIO declarations
cdef extern from "<OpenMS/FORMAT/ProteinIdentificationArrowIO.h>" namespace "OpenMS":

    cdef cppclass ProteinIdentificationArrowIO:
        # wrap-doc:
        #  Import and export ProteinIdentification data to/from Apache Arrow/Parquet format.
        #
        #  Provides static methods to export and import protein hits, protein groups, and search
        #  parameters to/from Parquet files. Enables full round-trip serialization.
        #
        #  EXPERIMENTAL: This API is experimental and may change in future versions.
        ProteinIdentificationArrowIO() except + nogil
        ProteinIdentificationArrowIO(const ProteinIdentificationArrowIO&) except + nogil

        bool exportProteinsToParquet(
            const libcpp_vector[ProteinIdentification]& protein_identifications,
            const String& filename,
            const ParquetWriteConfig& config
        ) except + nogil
        # wrap-doc:
        #  Export protein hits to Parquet file (one row per ProteinHit).
        #  Returns true on success, false on error.

        bool exportProteinGroupsToParquet(
            const libcpp_vector[ProteinIdentification]& protein_identifications,
            const String& filename,
            const ParquetWriteConfig& config
        ) except + nogil
        # wrap-doc:
        #  Export protein groups to Parquet file (one row per group).
        #  Returns true on success, false on error.

        bool exportSearchParamsToParquet(
            const libcpp_vector[ProteinIdentification]& protein_identifications,
            const String& filename,
            const ParquetWriteConfig& config
        ) except + nogil
        # wrap-doc:
        #  Export search parameters to Parquet file (one row per search run).
        #  Returns true on success, false on error.

        bool importFromParquet(
            const String& proteins_filename,
            const String& protein_groups_filename,
            const String& search_params_filename,
            libcpp_vector[ProteinIdentification]& protein_identifications
        ) except + nogil
        # wrap-doc:
        #  Import all three Parquet files and reconstruct ProteinIdentifications.
        #  Returns true on success, false on error.

        bool importSearchParamsFromParquet(
            const String& filename,
            libcpp_vector[ProteinIdentification]& protein_identifications
        ) except + nogil
        # wrap-doc:
        #  Import search parameters from Parquet file.
        #  Returns true on success, false on error.

        bool importProteinsFromParquet(
            const String& filename,
            libcpp_vector[ProteinIdentification]& protein_identifications
        ) except + nogil
        # wrap-doc:
        #  Import protein hits from Parquet file.
        #  Returns true on success, false on error.

        bool importProteinGroupsFromParquet(
            const String& filename,
            libcpp_vector[ProteinIdentification]& protein_identifications
        ) except + nogil
        # wrap-doc:
        #  Import protein groups from Parquet file.
        #  Returns true on success, false on error.
