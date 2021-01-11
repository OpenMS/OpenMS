from MSExperiment cimport *
from MSChromatogram cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

from MassTrace cimport *
from Feature cimport *
from FeatureMap cimport *
from MSChromatogram cimport *

from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>" namespace "OpenMS":

    cdef cppclass FeatureFindingMetabo(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #    ProgressLogger
        #    DefaultParamHandler
        #

        FeatureFindingMetabo()      nogil except +

        void run(libcpp_vector[Kernel_MassTrace] input_mtraces,
                 FeatureMap & output_featmap,
                 libcpp_vector[ libcpp_vector[ MSChromatogram ] ] & output_chromatograms
                 ) nogil except +

