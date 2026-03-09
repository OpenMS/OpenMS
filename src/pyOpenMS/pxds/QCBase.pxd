from libcpp cimport bool
from Types cimport *
from String cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/QC/QCBase.h>" namespace "OpenMS":

    cdef cppclass QCBase:
        QCBase() nogil except +
        String getName() nogil except +
        bool isRunnable(QCBase.Status) nogil except +

        cdef cppclass SpectraMap:
            SpectraMap() nogil except +
            void calculateMap(MSExperiment) nogil except +
            UInt64 at(String) nogil except +
            void clear() nogil except +
            bool empty() nogil except +
            Size size() nogil except +