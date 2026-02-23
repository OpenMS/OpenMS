from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport uint32_t, uint64_t, int64_t
from OpenMS cimport *
from QCBase cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from PeptideIdentification cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/QC/Ms2SpectrumStats.h>" namespace "OpenMS":

    cdef cppclass Ms2SpectrumStats(QCBase):
        # wrap-inherits:
        #    QCBase

        # Nested struct
        cppclass ScanEvent:
            ScanEvent(uint32_t sem, bool ms2) except +
            uint32_t scan_event_number
            bool ms2_presence

        Ms2SpectrumStats() except +
        Ms2SpectrumStats(Ms2SpectrumStats &) except +

        PeptideIdentificationList compute(const MSExperiment& exp, 
                                         FeatureMap& features, 
                                         const QCBase.SpectraMap& map_to_spectrum) except +

        String getName() const

        QCBase.Status requirements() const
