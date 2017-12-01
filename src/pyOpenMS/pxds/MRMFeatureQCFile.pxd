from String cimport *
from MRMFeatureQC cimport *

cdef extern from "<OpenMS/FORMAT/MRMFeatureQCFile.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureQCFile:

        MRMFeatureQCFile() nogil except +
        MRMFeatureQCFile(MRMFeatureQCFile &) nogil except +

        void load(String filename, MRMFeatureQC mrmfqc) nogil except +
        # void store(String filename, MRMFeatureQC mrmfqc) nogil except +

