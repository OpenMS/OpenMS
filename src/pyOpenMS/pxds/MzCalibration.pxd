from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport uint32_t, uint64_t, int64_t
from OpenMS cimport *
from QCBase cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/QC/MzCalibration.h>" namespace "OpenMS":

    cdef cppclass MzCalibration(QCBase):
        # wrap-inherits:
        #    QCBase

        MzCalibration() except +
        MzCalibration(MzCalibration &) except +

        void compute(FeatureMap& features, 
                     const MSExperiment& exp, 
                     const QCBase.SpectraMap& map_to_spectrum) except +

        String getName() except + const 

        QCBase.Status requirements() except + const
