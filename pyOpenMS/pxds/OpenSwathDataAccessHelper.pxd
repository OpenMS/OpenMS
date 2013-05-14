from MSSpectrum cimport *
from LightTargetedExperiment cimport LightTargetedExperiment
from TargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>" namespace "OpenMS":

    cdef cppclass OpenSwathDataAccessHelper:

        void convertTargetedExp(TargetedExperiment & transition_exp_, LightTargetedExperiment & transition_exp)

