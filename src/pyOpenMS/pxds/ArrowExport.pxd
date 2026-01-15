# Part of pyOpenMS. Copyright 2002-present, OpenMS Inc.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from Types cimport *
from MSExperiment cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from libcpp.string cimport string as libcpp_utf8_output_string
from libcpp cimport bool
from libc.stdint cimport int64_t, uint8_t

# OpenMS ArrowExport declarations (only available when WITH_PARQUET is enabled)
cdef extern from "<OpenMS/FORMAT/ArrowExport.h>" namespace "OpenMS":

    cdef cppclass ArrowSpectraExportConfig:
        # wrap-doc:
        #  Configuration for Arrow export of spectra data.
        #
        #  Allows filtering by MS level, RT range, m/z range, and column selection.
        #  Set ms_levels to filter specific MS levels (empty = all levels).
        #  Set min_rt/max_rt and min_mz/max_mz for range filtering (0 = no filter).
        #  Set columns to export specific columns (empty = all columns).
        ArrowSpectraExportConfig() except + nogil
        ArrowSpectraExportConfig(const ArrowSpectraExportConfig&) except + nogil  # copy constructor
        # ArrowExportFormat format  # wrap-ignore (enum class not supported)
        libcpp_vector[unsigned int] ms_levels
        double min_rt
        double max_rt
        double min_mz
        double max_mz
        libcpp_vector[libcpp_string] columns
        bool include_precursor_info
        bool include_ion_mobility

    cdef cppclass ArrowChromatogramExportConfig:
        # wrap-doc:
        #  Configuration for Arrow export of chromatogram data.
        #
        #  Allows filtering by RT range and column selection.
        #  Set min_rt/max_rt for range filtering (0 = no filter).
        #  Set columns to export specific columns (empty = all columns).
        ArrowChromatogramExportConfig() except + nogil
        ArrowChromatogramExportConfig(const ArrowChromatogramExportConfig&) except + nogil  # copy constructor
        # ArrowExportFormat format  # wrap-ignore (enum class not supported)
        double min_rt
        double max_rt
        libcpp_vector[libcpp_string] columns

    cdef cppclass ParquetWriteConfig:
        # wrap-doc:
        #  Configuration for Parquet file writing.
        #
        #  Controls compression, row group size, and other Parquet-specific settings.
        #  Default settings use ZSTD compression with 128MB row groups.
        ParquetWriteConfig() except + nogil
        ParquetWriteConfig(const ParquetWriteConfig&) except + nogil  # copy constructor
        # Compression compression  # wrap-ignore (enum class not supported, defaults to ZSTD)
        int compression_level
        int64_t row_group_size
        bool write_statistics
        int64_t data_page_size

    cdef cppclass ArrowExport:
        # wrap-doc:
        #  Export MSExperiment data to Apache Arrow format.
        #
        #  This class provides static methods to export MSExperiment spectra and
        #  chromatograms to Apache Arrow Tables and Parquet files. The export is
        #  efficient for large datasets as it uses Arrow's columnar memory format.
        #
        #  EXPERIMENTAL: This API is experimental and may change in future versions.
        #  The table schema, column names, and data types are subject to modification.
        ArrowExport() except + nogil
        ArrowExport(const ArrowExport&) except + nogil  # copy constructor
        # Column discovery - these are static methods but declared as instance methods for autowrap
        libcpp_vector[libcpp_utf8_output_string] getSpectraArrowColumnNames(
            const MSExperiment& exp,
            const ArrowSpectraExportConfig& config
        ) except + nogil
        # wrap-doc:
        #  Get available column names for spectra Arrow export.
        #  This is a lightweight operation that returns the list of columns that would
        #  be included in the export based on the configuration. Useful for discovering
        #  columns before calling to_arrow() or for column selection.

        libcpp_vector[libcpp_utf8_output_string] getChromatogramArrowColumnNames(
            const MSExperiment& exp,
            const ArrowChromatogramExportConfig& config
        ) except + nogil
        # wrap-doc:
        #  Get available column names for chromatogram Arrow export.
        #  This is a lightweight operation that returns the list of columns that would
        #  be included in the export based on the configuration.

        # Parquet export - these are static methods but declared as instance methods for autowrap
        bool exportSpectraToParquet(
            const MSExperiment& exp,
            const String& filename,
            const ArrowSpectraExportConfig& config,
            const ParquetWriteConfig& parquet_config
        ) except + nogil
        # wrap-doc:
        #  Export spectra to Parquet file. Parquet provides columnar storage with
        #  efficient compression (typically 3-5x for MS data). The operation iterates
        #  through all spectra matching the filter criteria, so runtime is O(n) where
        #  n is the number of filtered peaks.

        bool exportChromatogramsToParquet(
            const MSExperiment& exp,
            const String& filename,
            const ArrowChromatogramExportConfig& config,
            const ParquetWriteConfig& parquet_config
        ) except + nogil
        # wrap-doc:
        #  Export chromatograms to Parquet file. Similar to exportSpectraToParquet
        #  but for chromatogram data. Runtime is O(n) where n is the number of
        #  filtered data points.
