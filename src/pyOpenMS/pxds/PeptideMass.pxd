from libcpp.vector cimport vector
from libcpp.string cimport string
from OpenMS cimport *
from QCBase cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/QC/PeptideMass.h>" namespace "OpenMS":

    cdef cppclass PeptideMass(QCBase):
        # wrap-inherits:
        #    QCBase

        PeptideMass() except +
        PeptideMass(PeptideMass &) except +

        void compute(FeatureMap& features) except +

        String getName() except + const

        QCBase.Status requirements() except + const
