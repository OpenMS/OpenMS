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

cdef extern from "<OpenMS/FEATUREFINDER/FeatureFindingMetabo.h>" namespace "OpenMS":

    cdef cppclass FeatureFindingMetabo(ProgressLogger, DefaultParamHandler):
        # wrap-inherits:
        #   ProgressLogger
        #   DefaultParamHandler
        #

        FeatureFindingMetabo() except + nogil 
            # wrap-doc: 
            #  Method for the assembly of mass traces belonging to the same isotope
            #  pattern, i.e., that are compatible in retention times, mass-to-charge ratios,
            #  and isotope abundances

        FeatureFindingMetabo(FeatureFindingMetabo &) except + nogil  # compiler

        void run(libcpp_vector[Kernel_MassTrace] input_mtraces,
                 FeatureMap & output_featmap,
                 libcpp_vector[ libcpp_vector[ MSChromatogram ] ] & output_chromatograms
                 ) except + nogil 

