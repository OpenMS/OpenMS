from Types cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from FloatDataArray cimport *
from IMTypes cimport *
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "<OpenMS/IONMOBILITY/IMDataConverter.h>" namespace "OpenMS":
    
    cdef cppclass IMDataConverter:

        IMDataConverter() except + nogil  # compiler
        IMDataConverter(IMDataConverter &) except + nogil  # compiler

# COMMENT: wrap static methods
cdef extern from "<OpenMS/IONMOBILITY/IMDataConverter.h>" namespace "OpenMS::IMDataConverter":

        vector[pair[double, MSExperiment]] splitByFAIMSCV(MSExperiment exp) except + nogil  # wrap-attach:IMDataConverter

        MSExperiment reshapeIMFrameToMany(MSSpectrum im_frame) except + nogil  # wrap-attach:IMDataConverter

        MSExperiment reshapeIMFrameToSingle(MSExperiment exp) except + nogil  # wrap-attach:IMDataConverter

        void setIMUnit(FloatDataArray& fda, DriftTimeUnit unit) except + nogil  # wrap-attach:IMDataConverter

        bool getIMUnit(FloatDataArray& fda, DriftTimeUnit& unit) except + nogil  # wrap-attach:IMDataConverter

        double convertVSSCToCCS(double IM, double mz, int charge) except + nogil  # wrap-attach:IMDataConverter wrap-as:convertVSSCToCCSSingle

        void convertVSSCToCCS(MSExperiment& spectra) except + nogil  # wrap-attach:IMDataConverter

