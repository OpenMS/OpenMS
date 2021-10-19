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
        # wrap-doc:
        #   This class implements the database preprocessing needing for precursor ion selection

        PrecursorIonSelectionPreprocessing() nogil except +
        PrecursorIonSelectionPreprocessing(PrecursorIonSelectionPreprocessing &) nogil except +

        # TODO STL map
        # libcpp_map[ String, libcpp_vector[ double ] ]  getProtMasses() nogil except +
        libcpp_vector[ double ]  getMasses(String acc) nogil except +
        # libcpp_map[ String, libcpp_vector[ double ] ]  getProteinRTMap() nogil except +
        # libcpp_map[ String, libcpp_vector[ double ] ]  getProteinPTMap() nogil except +
        # libcpp_map[ String, libcpp_vector[ String ] ]  getProteinPeptideSequenceMap() nogil except +
        void dbPreprocessing(String db_path, bool save) nogil except +
            # wrap-doc:
                #   Calculates tryptic peptide masses of a given database and stores masses and peptide sequences
                #   -----
                #   :param db_path: Path to database file (fasta)
                #   :param save: Flag if preprocessing should be stored
                #   -----
                #   :raises:
                #     Exception: FileNotFound is thrown if the file could not be found
                #   :raises:
                #     Exception: UnableToCreateFile if preprocessing file can't be written

        void dbPreprocessing(String db_path, String rt_model_path, String dt_model_path, bool save) nogil except +
            # wrap-doc:
                #   Calculates tryptic peptide masses of a given database and stores masses and peptide sequences
                #   -----
                #   :param db_path: Path to database file (fasta)
                #   :param rt_model_path
                #   :param dt_model_path
                #   :param save: Flag if preprocessing should be stored
                #   -----
                #   :raises:
                #     Exception: FileNotFound is thrown if the file could not be found
                #   :raises:
                #     Exception: UnableToCreateFile if preprocessing file can't be written

        void loadPreprocessing() nogil except + # wrap-doc:Loads tryptic peptide masses of a given database
        double getWeight(double mass) nogil except + # wrap-doc:Returns the weighted frequency of a mass
        double getRT(String prot_id, Size peptide_index) nogil except + # wrap-doc:Returns the RT value
        double getPT(String prot_id, Size peptide_index) nogil except + # wrap-doc:Returns the PT value
        void setFixedModifications(StringList & modifications) nogil except +
        # libcpp_map[ char, libcpp_vector[ String ] ]  getFixedModifications() nogil except +
        void setGaussianParameters(double mu, double sigma) nogil except +
        double getGaussMu() nogil except + # wrap-doc:Returns the Gauss Mu value
        double getGaussSigma() nogil except + # wrap-doc:Returns the Gauss Sigma value
        double getRTProbability(String prot_id, Size peptide_index, Feature & feature) nogil except +
        double getRTProbability(double pred_rt, Feature & feature) nogil except +

