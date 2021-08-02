from String cimport *
from MRMFeatureQC cimport *

cdef extern from "<OpenMS/FORMAT/MRMFeatureQCFile.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureQCFile:
        # wrap-doc:
                #   File adapter for MRMFeatureQC files
                #   -----
                #   Loads and stores .csv or .tsv files describing an MRMFeatureQC

        MRMFeatureQCFile() nogil except +
        MRMFeatureQCFile(MRMFeatureQCFile &) nogil except + # compiler

        void load(const String& filename, MRMFeatureQC& mrmfqc, const bool is_component_group) nogil except +
            # wrap-doc:
                #   Loads an MRMFeatureQC file
                #   -----
                #   :param filename: The path to the input file
                #   :param mrmfqc: The output class which will contain the criteria
                #   :param is_component_group: True if the user intends to load ComponentGroupQCs data, false otherwise
                #   -----
                #   :raises:
                #     Exception: FileNotFound is thrown if the file could not be opened
                #   :raises:
                #     Exception: ParseError is thrown if an error occurs during parsing

        # void store(String filename, MRMFeatureQC mrmfqc) nogil except +
