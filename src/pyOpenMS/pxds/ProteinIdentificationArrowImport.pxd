# Part of pyOpenMS. Copyright 2002-present, OpenMS Inc.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from Types cimport *
from ProteinIdentification cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

# OpenMS ProteinIdentificationArrowImport declarations (only available when WITH_PARQUET is enabled)
cdef extern from "<OpenMS/FORMAT/ProteinIdentificationArrowImport.h>" namespace "OpenMS":

    cdef cppclass ProteinIdentificationArrowImport:
        # wrap-doc:
        #  Import ProteinIdentification data from Apache Arrow/Parquet format.
        #
        #  Provides static methods to import protein hits, protein groups, and search
        #  parameters from Parquet files. Reverse of ProteinIdentificationArrowExport.
        #
        #  EXPERIMENTAL: This API is experimental and may change in future versions.
        ProteinIdentificationArrowImport() except + nogil
        ProteinIdentificationArrowImport(const ProteinIdentificationArrowImport&) except + nogil

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
