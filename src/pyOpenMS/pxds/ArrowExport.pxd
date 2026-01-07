# Part of pyOpenMS. Copyright 2002-present, OpenMS Inc.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause

from Types cimport *
from MSExperiment cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from libcpp cimport bool
from libc.stdint cimport int64_t, uint8_t

# OpenMS ArrowExport declarations (only available when WITH_PARQUET is enabled)
cdef extern from "<OpenMS/FORMAT/ArrowExport.h>" namespace "OpenMS":

    cdef cppclass ArrowSpectraExportConfig:
        # wrap-doc:
        #  Configuration for Arrow export of spectra data.
        #
        #  Allows filtering by MS level, RT range, m/z range, and column selection.
        ArrowSpectraExportConfig() except + nogil
        ArrowSpectraExportConfig(ArrowSpectraExportConfig&) except + nogil  # copy constructor
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
        ArrowChromatogramExportConfig() except + nogil
        ArrowChromatogramExportConfig(ArrowChromatogramExportConfig&) except + nogil  # copy constructor
        # ArrowExportFormat format  # wrap-ignore (enum class not supported)
        double min_rt
        double max_rt
        libcpp_vector[libcpp_string] columns

    cdef cppclass ArrowExport:
        # wrap-doc:
        #  Export MSExperiment data to Apache Arrow format.
        #
        #  This class provides static methods to export MSExperiment spectra and
        #  chromatograms to Apache Arrow Tables and Parquet files.
        ArrowExport() except + nogil

        # Column discovery - these are static methods but declared as instance methods for autowrap
        libcpp_vector[libcpp_string] getSpectraArrowColumns(
            MSExperiment& exp,
            ArrowSpectraExportConfig& config
        ) except + nogil  # wrap-doc:Get available column names for spectra Arrow export

        libcpp_vector[libcpp_string] getChromatogramArrowColumns(
            MSExperiment& exp,
            ArrowChromatogramExportConfig& config
        ) except + nogil  # wrap-doc:Get available column names for chromatogram Arrow export

        # Parquet export - these are static methods but declared as instance methods for autowrap
        bool exportSpectraToParquet(
            MSExperiment& exp,
            const String& filename,
            ArrowSpectraExportConfig& config
        ) except + nogil  # wrap-doc:Export spectra to Parquet file

        bool exportChromatogramsToParquet(
            MSExperiment& exp,
            const String& filename,
            ArrowChromatogramExportConfig& config
        ) except + nogil  # wrap-doc:Export chromatograms to Parquet file
