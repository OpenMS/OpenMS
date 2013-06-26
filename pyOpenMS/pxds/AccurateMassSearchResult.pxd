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
        DoubleReal getAdductMass() nogil except +
        void setAdductMass(DoubleReal & m) nogil except +
        DoubleReal getQueryMass() nogil except +
        void setQueryMass(DoubleReal & m) nogil except +
        DoubleReal getFoundMass() nogil except +
        void setFoundMass(DoubleReal & m) nogil except +
        DoubleReal getCharge() nogil except +
        void setCharge(DoubleReal & ch) nogil except +
        DoubleReal getErrorPPM() nogil except +
        void setErrorPPM(DoubleReal & ppm) nogil except +
        DoubleReal getObservedRT() nogil except +
        void setObservedRT(DoubleReal & rt) nogil except +
        DoubleReal getObservedIntensity() nogil except +
        void setObservedIntensity(DoubleReal & intensity) nogil except +
        DoubleReal getMatchingIndex() nogil except +
        void setMatchingIndex(DoubleReal & idx) nogil except +
        String getFoundAdduct() nogil except +
        void setFoundAdduct(String & add) nogil except +
        String getFormulaString() nogil except +
        void setEmpiricalFormula(String & ep) nogil except +
        libcpp_vector[ String ] getMatchingHMDBids() nogil except +
        void setMatchingHMDBids(libcpp_vector[ String ] & match_ids) nogil except +
        DoubleReal getIsotopesSimScore() nogil except +
        void setIsotopesSimScore(DoubleReal & sim_score) nogil except +
        void outputResults() nogil except +

