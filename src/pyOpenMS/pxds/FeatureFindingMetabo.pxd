from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

from MassTrace cimport *
from Feature cimport *
from FeatureMap cimport *

from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>" namespace "OpenMS":

    cdef cppclass FeatureFindingMetabo(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #    ProgressLogger
        #    DefaultParamHandler
        #

        FeatureFindingMetabo()      nogil except +

        void run(libcpp_vector[Kernel_MassTrace] input,
                 FeatureMap[Feature] & result
                 ) nogil except +

