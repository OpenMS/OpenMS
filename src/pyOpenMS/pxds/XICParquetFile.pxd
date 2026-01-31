from libcpp.vector cimport vector
from libcpp cimport bool

from String cimport String
from StringList cimport StringList
from Types cimport Int64

ctypedef vector[double] DoubleVector

cdef extern from "<OpenMS/FORMAT/XICParquetFile.h>" namespace "OpenMS":
    cdef cppclass XICChromatogram:
        # wrap-ignore
        Int64 run_id
        String source_file
        Int64 ms_level

        bool has_precursor_id
        Int64 precursor_id
        bool has_transition_id
        Int64 transition_id
        String modified_sequence
        bool has_precursor_charge
        Int64 precursor_charge
        bool has_product_charge
        Int64 product_charge
        bool has_detecting_transition
        Int64 detecting_transition
        bool has_precursor_decoy
        Int64 precursor_decoy
        bool has_product_decoy
        Int64 product_decoy
        bool has_transition_ordinal
        Int64 transition_ordinal
        String transition_type
        String annotation

        DoubleVector rt  # wrap-ignore
        DoubleVector intensity  # wrap-ignore

    cdef cppclass XICParquetFile:
        XICParquetFile(String filename) except + nogil
            # wrap-doc:
            #  Reader for OpenSWATH chromatogram Parquet files (.xic).
        XICParquetFile(const StringList& filenames) except + nogil
            # wrap-doc:
            #  Reader for multiple OpenSWATH chromatogram Parquet files (.xic).

        const String& getFilename() const
        const StringList& getFilenames() const

        void load(vector[XICChromatogram]& output) const  # wrap-ignore
            # wrap-doc:
            #  Load all chromatograms from the file.

        void getChromatograms(vector[XICChromatogram]& output,
                              Int64 precursor_id,
                              Int64 transition_id,
                              String modified_sequence,
                              Int64 precursor_charge,
                              Int64 product_charge,
                              Int64 ms_level,
                              Int64 run_id,
                              String filter) const  # wrap-ignore
            # wrap-doc:
            #  Load chromatograms with optional filtering.
            #  
            #  :param precursor_id: Optional precursor id (-1 to ignore)
            #  :param transition_id: Optional transition id (-1 to ignore)
            #  :param modified_sequence: Optional modified sequence (empty to ignore)
            #  :param precursor_charge: Optional precursor charge (-1 to ignore)
            #  :param product_charge: Optional product charge (-1 to ignore)
            #  :param ms_level: Optional MS level (-1 to ignore)
            #  :param run_id: Optional run id (-1 to ignore)
            #  :param filter: Optional filter expression string
