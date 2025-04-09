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
        AccurateMassSearchResult() except + nogil 
        AccurateMassSearchResult(AccurateMassSearchResult &) except + nogil  # wrap-ignore
        double getObservedMZ() except + nogil 
        void setObservedMZ(double & m) except + nogil 
        double getCalculatedMZ() except + nogil 
        void setCalculatedMZ(double & m) except + nogil 
        double getQueryMass() except + nogil 
        void setQueryMass(double & m) except + nogil 
        double getFoundMass() except + nogil 
        void setFoundMass(double & m) except + nogil 
        double getCharge() except + nogil 
        void setCharge(double & ch) except + nogil 
        double getMZErrorPPM() except + nogil 
        void setMZErrorPPM(double & ppm) except + nogil 
        double getObservedRT() except + nogil 
        void setObservedRT(double & rt) except + nogil 
        double getObservedIntensity() except + nogil 
        void setObservedIntensity(double & intensity) except + nogil 
        double getMatchingIndex() except + nogil 
        void setMatchingIndex(double & idx) except + nogil 
        String getFoundAdduct() except + nogil 
        void setFoundAdduct(const String & add) except + nogil 
        String getFormulaString() except + nogil 
        void setEmpiricalFormula(const String & ep) except + nogil 
        libcpp_vector[ String ] getMatchingHMDBids() except + nogil 
        void setMatchingHMDBids(libcpp_vector[ String ] & match_ids) except + nogil 
        double getIsotopesSimScore() except + nogil 
        void setIsotopesSimScore(double & sim_score) except + nogil 
        libcpp_vector[double] getIndividualIntensities() except + nogil 
        void setIndividualIntensities(libcpp_vector[double]) except + nogil 
        Size getSourceFeatureIndex() except + nogil 
        void setSourceFeatureIndex(Size) except + nogil 
        # errors when returning references to vector
        libcpp_vector[ double] getMasstraceIntensities() except + nogil 
        void setMasstraceIntensities(libcpp_vector[ double ] & ) except + nogil 
