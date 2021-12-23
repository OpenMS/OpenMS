from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.string cimport string as libcpp_utf8_string
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/FORMAT/OSWFile.h>" namespace "OpenMS":
    
    cdef cppclass OSWFile "OpenMS::OSWFile":
        # wrap-doc:
            #   This class serves for reading in and writing OpenSWATH OSW files
            #   -----
            #   See OpenSwathOSWWriter for more functionality
            #   -----
            #   The reader and writer returns data in a format suitable for PercolatorAdapter.
            #   OSW files have a flexible data structure. They contain all peptide query
            #   parameters of TraML/PQP files with the detected and quantified features of
            #   OpenSwathWorkflow (feature, feature_ms1, feature_ms2 & feature_transition)
            #   -----
            #   The OSWFile reader extracts the feature information from the OSW file for
            #   each level (MS1, MS2 & transition) separately and generates Percolator input
            #   files. For each of the three Percolator reports, OSWFile writer adds a table
            #   (score_ms1, score_ms2, score_transition) with the respective confidence metrics.
            #   These tables can be mapped to the corresponding feature tables, are very similar
            #   to PyProphet results and can thus be used interchangeably

        OSWFile(const libcpp_utf8_string filename) nogil except +
        OSWFile(OSWFile &) nogil except +
        
        # Cannot wrap libcpp_ostream
        # void readToPIN(const libcpp_string & in_osw,
        #          const int osw_level,
        #          libcpp_ostream & pin_output,
        #          const double ipf_max_peakgroup_pep,
        #          const double ipf_max_transition_isotope_overlap,
        #          const double ipf_min_transition_sn) nogil except +

        # NESTED STL 
        # void writeFromPercolator(const libcpp_string & in_osw, const int osw_level, const libcpp_map[ libcpp_string, PercolatorFeature ] & features) nogil except +

