from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from AASequence cimport *
from FASTAFile cimport *
from FeatureMap cimport *
from StringList cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>" namespace "OpenMS":
    
    cdef cppclass PrecursorIonSelectionPreprocessing(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PrecursorIonSelectionPreprocessing() nogil except +
        PrecursorIonSelectionPreprocessing(PrecursorIonSelectionPreprocessing) nogil except +

        # TODO STL map
        # libcpp_map[ String, libcpp_vector[ double ] ]  getProtMasses() nogil except +
        libcpp_vector[ double ]  getMasses(String acc) nogil except +
        # libcpp_map[ String, libcpp_vector[ double ] ]  getProteinRTMap() nogil except +
        # libcpp_map[ String, libcpp_vector[ double ] ]  getProteinPTMap() nogil except +
        # libcpp_map[ String, libcpp_vector[ String ] ]  getProteinPeptideSequenceMap() nogil except +
        void dbPreprocessing(String db_path, bool save) nogil except +
        void dbPreprocessing(String db_path, String rt_model_path, String dt_model_path, bool save) nogil except +
        void loadPreprocessing() nogil except +
        double getWeight(double mass) nogil except +
        double getRT(String prot_id, Size peptide_index) nogil except +
        double getPT(String prot_id, Size peptide_index) nogil except +
        void setFixedModifications(StringList & modifications) nogil except +
        # libcpp_map[ char, libcpp_vector[ String ] ]  getFixedModifications() nogil except +
        void setGaussianParameters(double mu, double sigma) nogil except +
        double getGaussMu() nogil except +
        double getGaussSigma() nogil except +
        double getRTProbability(String prot_id, Size peptide_index, Feature & feature) nogil except +
        double getRTProbability(double pred_rt, Feature & feature) nogil except +

