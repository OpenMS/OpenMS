from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/QC/MzCalibration.h>" namespace "OpenMS":
    
    cdef cppclass MzCalibration(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  QC metric calculating (un)calibrated m/z error
        #  
        #  The metric sets m/z-values of the original experiment and the calculated 
        #  reference m/z-values, uncalibrated m/z error (ppm) and calibrated m/z error (ppm) 
        #  as metavalues of all PeptideIdentifications in a FeatureMap.

        MzCalibration() except + nogil 
        MzCalibration(MzCalibration &) except + nogil 

        void compute(FeatureMap & features, MSExperiment & exp, SpectraMap & map_to_spectrum) except + nogil 
            # wrap-doc:
            #  Writes results as meta values to the PeptideIdentification of the given FeatureMap
            #  
            #  :param features: FeatureMap with m/z-values of PeptideIdentification after calibration
            #  :param exp: PeakMap of the original experiment
            #  :param map_to_spectrum: Map to find index of spectrum given by meta value at PepID

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric
