# distutils: language = c++
# cython: language_level=3
"""
Zero-copy Arrow export from OpenMS MSExperiment.

This module provides zero-copy export of MS data to Apache Arrow format
using the Arrow C Data Interface. It is only available when OpenMS is
built with WITH_PARQUET=ON.

Usage:
    from pyopenms._arrow_zerocopy import spectra_to_arrow, chromatograms_to_arrow

    # Zero-copy export
    table = spectra_to_arrow(exp, format='long', ms_levels=[1, 2])
"""

from libc.stdlib cimport malloc, free
from libc.string cimport memset
from libc.stdint cimport int64_t, uint8_t, uintptr_t
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string

# Import MSExperiment
from MSExperiment cimport MSExperiment as _MSExperiment_cls
from MSExperiment cimport _MSExperiment

# Arrow C Data Interface structs
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


# OpenMS ArrowExport declarations
cdef extern from "<OpenMS/FORMAT/ArrowExport.h>" namespace "OpenMS":

    cdef enum ArrowExportFormat:
        Long "OpenMS::ArrowExportFormat::Long"
        SemiWide "OpenMS::ArrowExportFormat::SemiWide"

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

    bool exportSpectraToArrowCDataInterface(
        _MSExperiment& exp,
        ArrowSpectraExportConfig& config,
        ArrowSchema* out_schema,
        ArrowArray* out_array
    ) except + nogil

    bool exportChromatogramsToArrowCDataInterface(
        _MSExperiment& exp,
        ArrowChromatogramExportConfig& config,
        ArrowSchema* out_schema,
        ArrowArray* out_array
    ) except + nogil


def spectra_to_arrow(exp, format='long', ms_levels=None,
                     min_rt=0.0, max_rt=0.0, min_mz=0.0, max_mz=0.0,
                     columns=None, include_precursor_info=True,
                     include_ion_mobility=True):
    """
    Export spectra to Arrow Table using zero-copy C Data Interface.

    Parameters
    ----------
    exp : MSExperiment
        The experiment to export.
    format : str
        'long' (one row per peak) or 'semi_wide' (one row per spectrum).
    ms_levels : list of int, optional
        MS levels to include. None means all levels.
    min_rt, max_rt : float
        RT range filter. 0 means no filter.
    min_mz, max_mz : float
        m/z range filter. 0 means no filter.
    columns : list of str, optional
        Columns to include. None means all columns.
    include_precursor_info : bool
        Include precursor columns.
    include_ion_mobility : bool
        Include ion mobility column.

    Returns
    -------
    pyarrow.Table
        Arrow table with spectrum data (zero-copy from C++).
    """
    import pyarrow as pa

    # Get the C++ MSExperiment pointer
    cdef _MSExperiment* cpp_exp = (<_MSExperiment_cls?>exp).inst.get()

    # Configure export
    cdef ArrowSpectraExportConfig config
    config.format = ArrowExportFormat.Long if format == 'long' else ArrowExportFormat.SemiWide
    config.min_rt = min_rt
    config.max_rt = max_rt
    config.min_mz = min_mz
    config.max_mz = max_mz
    config.include_precursor_info = include_precursor_info
    config.include_ion_mobility = include_ion_mobility

    if ms_levels is not None:
        for level in ms_levels:
            config.ms_levels.push_back(<unsigned int>level)

    if columns is not None:
        for col in columns:
            config.columns.push_back(col.encode('utf-8'))

    # Allocate ArrowSchema and ArrowArray
    cdef ArrowSchema* schema = <ArrowSchema*>malloc(sizeof(ArrowSchema))
    cdef ArrowArray* array = <ArrowArray*>malloc(sizeof(ArrowArray))

    if schema == NULL or array == NULL:
        if schema != NULL:
            free(schema)
        if array != NULL:
            free(array)
        raise MemoryError("Failed to allocate Arrow structs")

    # Initialize to zero
    memset(schema, 0, sizeof(ArrowSchema))
    memset(array, 0, sizeof(ArrowArray))

    # Call C++ export
    cdef bint success
    with nogil:
        success = exportSpectraToArrowCDataInterface(cpp_exp[0], config, schema, array)

    if not success:
        free(schema)
        free(array)
        raise RuntimeError("Failed to export spectra to Arrow format")

    # Import into PyArrow using the C Data Interface
    # _import_from_c expects integer addresses (uintptr_t)
    # PyArrow takes ownership and will call the release callbacks
    try:
        batch = pa.RecordBatch._import_from_c(
            <uintptr_t>array,
            <uintptr_t>schema
        )
    except Exception as e:
        # On failure, we still own the memory - clean up
        if schema.release != NULL:
            schema.release(schema)
        if array.release != NULL:
            array.release(array)
        free(schema)
        free(array)
        raise RuntimeError(f"Failed to import Arrow data: {e}")

    # PyArrow now owns the data - free our allocations but not the Arrow data
    # (Arrow's release callbacks will be called by PyArrow when the table is garbage collected)
    free(schema)
    free(array)

    # Convert to Table
    return pa.Table.from_batches([batch])


def chromatograms_to_arrow(exp, format='long', min_rt=0.0, max_rt=0.0, columns=None):
    """
    Export chromatograms to Arrow Table using zero-copy C Data Interface.

    Parameters
    ----------
    exp : MSExperiment
        The experiment to export.
    format : str
        'long' (one row per point) or 'semi_wide' (one row per chromatogram).
    min_rt, max_rt : float
        RT range filter. 0 means no filter.
    columns : list of str, optional
        Columns to include. None means all columns.

    Returns
    -------
    pyarrow.Table
        Arrow table with chromatogram data (zero-copy from C++).
    """
    import pyarrow as pa

    # Get the C++ MSExperiment pointer
    cdef _MSExperiment* cpp_exp = (<_MSExperiment_cls?>exp).inst.get()

    # Configure export
    cdef ArrowChromatogramExportConfig config
    config.format = ArrowExportFormat.Long if format == 'long' else ArrowExportFormat.SemiWide
    config.min_rt = min_rt
    config.max_rt = max_rt

    if columns is not None:
        for col in columns:
            config.columns.push_back(col.encode('utf-8'))

    # Allocate ArrowSchema and ArrowArray
    cdef ArrowSchema* schema = <ArrowSchema*>malloc(sizeof(ArrowSchema))
    cdef ArrowArray* array = <ArrowArray*>malloc(sizeof(ArrowArray))

    if schema == NULL or array == NULL:
        if schema != NULL:
            free(schema)
        if array != NULL:
            free(array)
        raise MemoryError("Failed to allocate Arrow structs")

    # Initialize to zero
    memset(schema, 0, sizeof(ArrowSchema))
    memset(array, 0, sizeof(ArrowArray))

    # Call C++ export
    cdef bint success
    with nogil:
        success = exportChromatogramsToArrowCDataInterface(cpp_exp[0], config, schema, array)

    if not success:
        free(schema)
        free(array)
        raise RuntimeError("Failed to export chromatograms to Arrow format")

    # Import into PyArrow using the C Data Interface
    # _import_from_c expects integer addresses (uintptr_t)
    try:
        batch = pa.RecordBatch._import_from_c(
            <uintptr_t>array,
            <uintptr_t>schema
        )
    except Exception as e:
        # On failure, we still own the memory - clean up
        if schema.release != NULL:
            schema.release(schema)
        if array.release != NULL:
            array.release(array)
        free(schema)
        free(array)
        raise RuntimeError(f"Failed to import Arrow data: {e}")

    # PyArrow now owns the data - free our allocations
    free(schema)
    free(array)

    # Convert to Table
    return pa.Table.from_batches([batch])
