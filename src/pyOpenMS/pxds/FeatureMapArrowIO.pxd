# Part of pyOpenMS. Copyright 2002-present, OpenMS Inc.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from Types cimport *
from FeatureMap cimport *
from MSExperimentArrowExport cimport *
from libcpp cimport bool

# OpenMS FeatureMapArrowIO declarations
cdef extern from "<OpenMS/FORMAT/FeatureMapArrowIO.h>" namespace "OpenMS":

    cdef cppclass FeatureMapArrowIO:
        # wrap-doc:
        #  Import and export FeatureMap data to/from Apache Arrow/Parquet format.
        #
        #  Stores data as a directory of Parquet files (features, PSMs, proteins,
        #  protein groups, search parameters). Enables full round-trip serialization
        #  as a replacement for FeatureXML.
        #
        #  EXPERIMENTAL: This API is experimental and may change in future versions.
        FeatureMapArrowIO() except + nogil
        FeatureMapArrowIO(const FeatureMapArrowIO&) except + nogil

        bool exportToParquet(
            const FeatureMap& feature_map,
            const String& directory,
            const ParquetWriteConfig& config
        ) except + nogil
        # wrap-doc:
        #  Export FeatureMap to a directory of Parquet files.
        #  Creates features.parquet, psms.parquet, proteins.parquet,
        #  protein_groups.parquet, and search_params.parquet.
        #  Returns true on success, false on error.

        bool importFromParquet(
            const String& directory,
            FeatureMap& feature_map
        ) except + nogil
        # wrap-doc:
        #  Import FeatureMap from a directory of Parquet files.
        #  Reads all five Parquet files and reconstructs the complete FeatureMap.
        #  Returns true on success, false on error.
