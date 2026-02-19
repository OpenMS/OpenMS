# Part of pyOpenMS. Copyright 2002-present, OpenMS Inc.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from Types cimport *
from ProteinIdentification cimport *
from MSExperimentArrowExport cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

# OpenMS ProteinIdentificationArrowExport declarations (only available when WITH_PARQUET is enabled)
cdef extern from "<OpenMS/FORMAT/ProteinIdentificationArrowExport.h>" namespace "OpenMS":

    cdef cppclass ProteinIdentificationArrowExport:
        # wrap-doc:
        #  Export ProteinIdentification data to Apache Arrow/Parquet format.
        #
        #  Provides static methods to export protein hits, protein groups, and search
        #  parameters to Parquet files. Complements QPXFile PSM export for full idXML replacement.
        #
        #  EXPERIMENTAL: This API is experimental and may change in future versions.
        ProteinIdentificationArrowExport() except + nogil
        ProteinIdentificationArrowExport(const ProteinIdentificationArrowExport&) except + nogil

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
