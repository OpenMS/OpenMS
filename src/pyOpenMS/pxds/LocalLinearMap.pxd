from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from Matrix cimport *

cdef extern from "<OpenMS/ANALYSIS/PIP/LocalLinearMap.h>" namespace "OpenMS":
    
    cdef cppclass LocalLinearMap "OpenMS::LocalLinearMap":
        # wrap-doc:
            #   Trained Local Linear Map (LLM) model for peak intensity prediction
            #   -----
            #   This class offers a model for predictions of peptide peak heights
            #   (referred to as intensities) by a Local Linear Map (LLM) model and
            #   is the basis of PeakIntensityPredictor
            #   -----
            #   A general introduction to the Peak Intensity Predictor (PIP)
            #   can be found in the PIP Tutorial
            #   -----
            #   The model trained needs two files for storing the position of the
            #   codebook vectors and the linear mappings (codebooks.data, linearMapping.data)
            #   This is the default model used by PeakIntensityPredictor

        LocalLinearMap() nogil except +
        # private
        LocalLinearMap(LocalLinearMap &) nogil except + # wrap-ignore
        LLMParam getLLMParam() nogil except + # wrap-doc:Returns parameters of the LocalLinearMap model
        Matrix[ double ]  getCodebooks() nogil except + # wrap-doc:Returns position of the codebook vectors (18-dim)
        Matrix[ double ]  getMatrixA() nogil except + # wrap-doc:Returns linear mappings of the codebooks
        libcpp_vector[ double ]  getVectorWout() nogil except + # wrap-doc:Returns linear bias

        # TODO STL attributes unsigned int 
        # Matrix[ UInt ]  getCord() nogil except +
        void normalizeVector(libcpp_vector[ double ] & aaIndexVariables) nogil except + # wrap-doc:Calculates and returns the normalized amino acid index variables from string representation of peptide
        # libcpp_vector[ double ] neigh(Matrix[ unsigned int ] & cord, Size win, double radius) nogil except +


cdef extern from "<OpenMS/ANALYSIS/PIP/LocalLinearMap.h>" namespace "OpenMS::LocalLinearMap":
    
    cdef cppclass LLMParam "OpenMS::LocalLinearMap::LLMParam":
        LLMParam() nogil except +
        LLMParam(LLMParam) nogil except + #wrap-ignore
        UInt xdim
        UInt ydim
        double radius
