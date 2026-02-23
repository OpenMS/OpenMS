from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map as cppmap
from libc.stdint cimport uint32_t, uint64_t, int64_t
from OpenMS cimport *
from QCBase cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/QC/SpectrumCount.h>" namespace "OpenMS":

    cdef cppclass SpectrumCount(QCBase):
        # wrap-inherits:
        #    QCBase

        SpectrumCount() except +
        SpectrumCount(SpectrumCount &) except +

        cppmap[uint32_t, uint32_t] compute(const MSExperiment& exp) except +

        String getName() const

        QCBase.Status requirements() const
