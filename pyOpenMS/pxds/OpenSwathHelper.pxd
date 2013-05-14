from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from ProgressLogger cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>" namespace "OpenMS":

    cdef cppclass OpenSwathHelper:

        bool checkSwathMapAndSelectTransitions(MSExperiment[Peak1D, ChromatogramPeak] & exp, 
                                               TargetedExperiment & targeted_exp,
                                               TargetedExperiment & transition_exp_used,
                                               double min_upper_edge_dist)



