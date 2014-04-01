from Types cimport *
from MassTrace cimport *
from Feature cimport *
from FeatureMap cimport *
from MzTab cimport *
from DefaultParamHandler cimport *
from String cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>" namespace "OpenMS":
    
    cdef cppclass AccurateMassSearchResult "OpenMS::AccurateMassSearchResult":
        AccurateMassSearchResult() nogil except +
        AccurateMassSearchResult(AccurateMassSearchResult) nogil except +
        double getAdductMass() nogil except +
        void setAdductMass(double & m) nogil except +
        double getQueryMass() nogil except +
        void setQueryMass(double & m) nogil except +
        double getFoundMass() nogil except +
        void setFoundMass(double & m) nogil except +
        double getCharge() nogil except +
        void setCharge(double & ch) nogil except +
        double getErrorPPM() nogil except +
        void setErrorPPM(double & ppm) nogil except +
        double getObservedRT() nogil except +
        void setObservedRT(double & rt) nogil except +
        double getObservedIntensity() nogil except +
        void setObservedIntensity(double & intensity) nogil except +
        double getMatchingIndex() nogil except +
        void setMatchingIndex(double & idx) nogil except +
        String getFoundAdduct() nogil except +
        void setFoundAdduct(String & add) nogil except +
        String getFormulaString() nogil except +
        void setEmpiricalFormula(String & ep) nogil except +
        libcpp_vector[ String ] getMatchingHMDBids() nogil except +
        void setMatchingHMDBids(libcpp_vector[ String ] & match_ids) nogil except +
        double getIsotopesSimScore() nogil except +
        void setIsotopesSimScore(double & sim_score) nogil except +
        void outputResults() nogil except +

