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

# Arrow C Data Interface structs (from arrow/c/abi.h)
cdef extern from "<arrow/c/abi.h>" nogil:
    cdef struct ArrowSchema:
        const char* format
        const char* name
        const char* metadata
        int64_t flags
        int64_t n_children
        ArrowSchema** children
        ArrowSchema* dictionary
        void (*release)(ArrowSchema*)
        void* private_data

    cdef struct ArrowArray:
        int64_t length
        int64_t null_count
        int64_t offset
        int64_t n_buffers
        int64_t n_children
        const void** buffers
        ArrowArray** children
        ArrowArray* dictionary
        void (*release)(ArrowArray*)
        void* private_data


# OpenMS ArrowExport declarations (only available when WITH_PARQUET is enabled)
cdef extern from "<OpenMS/FORMAT/ArrowExport.h>" namespace "OpenMS":

    cdef cppclass ArrowExportFormat:
        pass

    ArrowExportFormat ArrowExportFormat_Long "OpenMS::ArrowExportFormat::Long"
    ArrowExportFormat ArrowExportFormat_SemiWide "OpenMS::ArrowExportFormat::SemiWide"

    cdef cppclass ArrowSpectraExportConfig:
        ArrowSpectraExportConfig() except + nogil
        ArrowExportFormat format
        libcpp_vector[unsigned int] ms_levels
        double min_rt
        double max_rt
        double min_mz
        double max_mz
        libcpp_vector[libcpp_string] columns
        bool include_precursor_info
        bool include_ion_mobility

    cdef cppclass ArrowChromatogramExportConfig:
        ArrowChromatogramExportConfig() except + nogil
        ArrowExportFormat format
        double min_rt
        double max_rt
        libcpp_vector[libcpp_string] columns

    # Column discovery functions
    libcpp_vector[libcpp_string] getSpectraArrowColumns(
        _MSExperiment& exp,
        ArrowSpectraExportConfig& config
    ) except + nogil  # wrap-doc:Get available column names for spectra Arrow export

    libcpp_vector[libcpp_string] getChromatogramArrowColumns(
        _MSExperiment& exp,
        ArrowChromatogramExportConfig& config
    ) except + nogil  # wrap-doc:Get available column names for chromatogram Arrow export

    # C Data Interface export functions (zero-copy to Python)
    bool exportSpectraToArrowCDataInterface(
        _MSExperiment& exp,
        ArrowSpectraExportConfig& config,
        ArrowSchema* out_schema,
        ArrowArray* out_array
    ) except + nogil  # wrap-doc:Export spectra to Arrow via C Data Interface

    bool exportChromatogramsToArrowCDataInterface(
        _MSExperiment& exp,
        ArrowChromatogramExportConfig& config,
        ArrowSchema* out_schema,
        ArrowArray* out_array
    ) except + nogil  # wrap-doc:Export chromatograms to Arrow via C Data Interface
