from libcpp cimport bool
from libcpp.vector cimport vector
from Types cimport *
from String cimport *
from MSExperiment cimport *
from MSChromatogram cimport *
from QCBase cimport *

cdef extern from "<OpenMS/QC/TIC.h>" namespace "OpenMS":

    cdef cppclass TIC(QCBase):
        TIC() nogil except +

        cdef cppclass Result:
            Result() nogil except +
            vector[UInt] intensities
            vector[float] relative_intensities
            vector[float] retention_times
            UInt area
            UInt fall
            UInt jump
            bool operator==(Result) nogil except +

        Result compute(MSExperiment, float, UInt) nogil except +
        String getName() nogil except +
        vector[MSChromatogram] getResults() nogil except +