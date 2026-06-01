# Part of pyOpenMS. Copyright 2002-present, OpenMS Inc.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from Types cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from MSExperimentArrowExport cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

# OpenMS QPXFile declarations
cdef extern from "<OpenMS/FORMAT/QPXFile.h>" namespace "OpenMS":

    cdef cppclass QPXFile:
        # wrap-doc:
        #  Export PSM (Peptide Spectrum Match) data to Apache Arrow format following QPX PSM schema.
        #
        #  This class provides static methods to export PeptideIdentification/ProteinIdentification
        #  data to Apache Arrow Tables and Parquet files. The schema follows the QPX (Quantitative
        #  Proteomics Exchange) PSM format.
        #
        #  EXPERIMENTAL: This API is experimental and may change in future versions.
        QPXFile() except + nogil
        QPXFile(const QPXFile&) except + nogil  # copy constructor

        bool exportToParquet(
            const libcpp_vector[ProteinIdentification]& protein_identifications,
            # PeptideIdentificationList not wrapped; use vector
            const libcpp_vector[PeptideIdentification]& peptide_identifications,
            const String& filename,
            bool export_all_psms,
            const ParquetWriteConfig& config
        ) except + nogil
        # wrap-doc:
        #  Export PSM data to Parquet file following the QPX PSM schema.
        #  Each PeptideHit becomes one row with identification, score, and metadata.
        #  Returns true on success, false on error.
