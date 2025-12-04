from Types cimport *
from DefaultParamHandler cimport *
from FeatureMap cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>" namespace "OpenMS":

    cdef cppclass Biosaur2Algorithm(DefaultParamHandler):
        # wrap-doc:
        #  C++ implementation of the Biosaur2 feature detection workflow.

        Biosaur2Algorithm() except + nogil  # wrap-doc:Instantiate the Biosaur2 feature finding algorithm
        Biosaur2Algorithm(Biosaur2Algorithm&) except + nogil  # wrap-doc:Copy constructor


        void setMSData(const MSExperiment& ms_data) except + nogil  # wrap-doc:Set the MS data used for feature detection (copy version)
        # wrap-ignore:
        #   void setMSData(MSExperiment&& ms_data) except + nogil
        
        MSExperiment& getMSData() except + nogil  # wrap-doc:Get non-const reference to MS data
        # wrap-as:getMSDataConst
        const MSExperiment& getMSData() except + nogil  # wrap-doc:Get const reference to MS data

        void run(FeatureMap& feature_map) except + nogil  # wrap-doc:Run the algorithm storing only the resulting features

