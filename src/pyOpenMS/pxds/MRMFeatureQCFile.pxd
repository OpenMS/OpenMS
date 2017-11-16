from String cimport *
from MRMFeatureQC import *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQCFile.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureQCFile:

        MRMFeatureQCFile() nogil except +
        MRMFeatureQCFile(MRMFeatureQCFile &) nogil except +

        void load(String filename, MRMFeatureQC mrmfqc) nogil except +
        void store(String filename, MRMFeatureQC mrmfqc);