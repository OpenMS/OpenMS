# Part of pyOpenMS. Copyright 2002-present, OpenMS Inc.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from Types cimport *
from ConsensusMap cimport *
from MSExperimentArrowExport cimport *
from libcpp cimport bool

# OpenMS ConsensusMapArrowIO declarations
cdef extern from "<OpenMS/FORMAT/ConsensusMapArrowIO.h>" namespace "OpenMS":

    cdef cppclass ConsensusMapArrowIO:
        # wrap-doc:
        #  Import and export ConsensusMap data to/from Apache Arrow/Parquet format.
        #
        #  Stores data as a directory of Parquet files (consensus features, PSMs,
        #  proteins, protein groups, search parameters). Enables full round-trip
        #  serialization as a replacement for ConsensusXML.
        #
        #  EXPERIMENTAL: This API is experimental and may change in future versions.
        ConsensusMapArrowIO() except + nogil
        ConsensusMapArrowIO(const ConsensusMapArrowIO&) except + nogil

        bool exportToParquet(
            const ConsensusMap& cmap,
            const String& directory,
            const ParquetWriteConfig& config
        ) except + nogil
        # wrap-doc:
        #  Export ConsensusMap to a directory of Parquet files.
        #  Creates consensus_features.parquet, psms.parquet, proteins.parquet,
        #  protein_groups.parquet, and search_params.parquet.
        #  Returns true on success, false on error.

        bool importFromParquet(
            const String& directory,
            ConsensusMap& cmap
        ) except + nogil
        # wrap-doc:
        #  Import ConsensusMap from a directory of Parquet files.
        #  Reads all five Parquet files and reconstructs the complete ConsensusMap.
        #  Returns true on success, false on error.
