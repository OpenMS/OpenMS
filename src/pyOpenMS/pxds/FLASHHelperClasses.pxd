from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from IsotopeDistribution cimport *
from Types cimport *



cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>" namespace "OpenMS":

    cdef cppclass MassFeature_FDHS "OpenMS::FLASHHelperClasses::MassFeature":
    
        # wrap-inherits:

        # default constructor
        MassFeature_FDHS() except + nogil
        # copy constructor
        MassFeature_FDHS(MassFeature_FDHS &) except + nogil


        bool operator<(MassFeature_FDHS& a) except + nogil
        bool operator>(MassFeature_FDHS& a) except + nogil
        bool operator==(MassFeature_FDHS& other) except + nogil

    cdef cppclass Tag "OpenMS::FLASHHelperClasses::Tag":

        # wrap-inherits:

        # default constructor
        Tag(String seq, double n_mass, double c_mass, libcpp_vector[double] & mzs, libcpp_vector[int]& scores, int scan) except + nogil
        # copy constructor
        Tag(Tag &) except + nogil

        String getSequence() except + nogil
        libcpp_vector[double] getMzs() except + nogil
        double getNtermMass() except + nogil
        double getCtermMass() except + nogil
        int getScore() except + nogil
        int getScore(int pos) except + nogil

    cdef cppclass PrecalAveragine "OpenMS::FLASHHelperClasses::PrecalculatedAveragine":
            
        # wrap-inherits:

        # default constructor
        PrecalAveragine() except + nogil
        # constructor
        PrecalAveragine(double min_mass, double max_mass, double delta, CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine) except + nogil
        # copy constructor
        PrecalAveragine(PrecalAveragine &) except + nogil
        
        IsotopeDistribution get(double mass) except + nogil
        size_t getMaxIsotopeIndex() except + nogil
        void setMaxIsotopeIndex(int index) except + nogil
        Size getLeftCountFromApex(double mass) except + nogil
        Size getRightCountFromApex(double mass) except + nogil
        Size getApexIndex(double mass) except + nogil
        Size getLastIndex(double mass) except + nogil
        double getAverageMassDelta(double mass) except + nogil
        double getMostAbundantMassDelta(double mass) except + nogil
        double getSNRMultiplicationFactor(double mass) except + nogil

    cdef cppclass IsobaricQuantities "OpenMS::FLASHHelperClasses::IsobaricQuantities":
       # wrap-inherits:

        # default constructor
        IsobaricQuantities() except + nogil
        # copy constructor
        IsobaricQuantities(IsobaricQuantities &) except + nogil

    cdef cppclass LogMzPeak "OpenMS::FLASHHelperClasses::LogMzPeak":
       # wrap-inherits:

        # default constructor
        LogMzPeak() except + nogil
        # copy constructor
        LogMzPeak(LogMzPeak &) except + nogil

        double getUnchargedMass() except + nogil

        bool operator<(LogMzPeak& a) except + nogil
        bool operator>(LogMzPeak& a) except + nogil
        bool operator==(LogMzPeak& other) except + nogil