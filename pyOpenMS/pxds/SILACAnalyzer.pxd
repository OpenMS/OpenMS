from ProgressLogger cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from String cimport *

from PeakWidthEstimator cimport *
from MSQuantifications cimport *

from ConsensusMap cimport *
from FeatureMap cimport *
from Feature cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SILACAnalyzer.h>" namespace "OpenMS":

    cdef cppclass SILACAnalyzer(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        SILACAnalyzer()                 nogil except +
        SILACAnalyzer(SILACAnalyzer)    nogil except + # wrap-ignore

        void initialize(String selected_labels_, UInt charge_min_, UInt charge_max_,
                        Int missed_cleavages_, UInt isotopes_per_peptide_min_, UInt isotopes_per_peptide_max_,
                        DoubleReal rt_threshold_, DoubleReal rt_min_, DoubleReal intensity_cutoff_,
                        DoubleReal intensity_correlation_, DoubleReal model_deviation_, 
                        bool allow_missing_peaks_, 
                        libcpp_map[String, DoubleReal] label_identifiers) nogil except + # wrap-ignore

        void run_all(MSExperiment[Peak1D, ChromatogramPeak] & exp, ConsensusMap & cmap) nogil except +

        void writeConsensus(String & filename, ConsensusMap & cmap) nogil except +

        # void calculateLabelsAndMassShifts(libcpp_map[ String, DoubleReal ] label_identifiers)
        Result estimatePeakWidth(MSExperiment[ Peak1D, ChromatogramPeak ] & exp)
        # TODO add SILACPattern, Clustering
        # NAMESPACE # void filterData(MSExperiment[ Peak1D, ChromatogramPeak ] & exp, PeakWidthEstimator::Result & peak_width, vector[ vector[ SILACPattern ] ] & data)
        # TODO nested STL
        # libcpp_vector[ libcpp_vector[ String ] ]  getSILAClabels()
        # libcpp_vector[ libcpp_vector[ double ] ]  getMassShifts()
        # void generateClusterConsensusByCluster(ConsensusMap & , Clustering & )
        # void generateClusterConsensusByPattern(ConsensusMap & , Clustering & , UInt & cluster_id)
        # void generateClusterDebug(std::ostream & out, Clustering & clustering, UInt & cluster_id)
        # void generateFilterConsensusByPattern(ConsensusMap & , libcpp_vector[ SILACPattern ] & )
        # void generateClusterFeatureByCluster(FeatureMap[Feature] & , Clustering & )
        # void readFilterConsensusByPattern(ConsensusMap & , vector[ vector[ SILACPattern ] ] & )
        void readConsensus(String & filename, ConsensusMap & in_)
        void writeMzQuantML(String & filename, MSQuantifications & msq)
        void writeFeatures(String & filename, FeatureMap[Feature] & out)
        String  selectColor(UInt nr)
        
