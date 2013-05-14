from String cimport *
from MSExperiment cimport *
from FeatureMap cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureFinderScoring:

        MRMFeatureFinderScoring() nogil except +

        void pickExperiment(MSExperiment[Peak1D, ChromatogramPeak] & chromatograms,
                            FeatureMap[Feature]& output,
                            TargetedExperiment& transition_exp_,
                            TransformationDescription trafo,
                            MSExperiment[Peak1D, ChromatogramPeak] & swath_map)


