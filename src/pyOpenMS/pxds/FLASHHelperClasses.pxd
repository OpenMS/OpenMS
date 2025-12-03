from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from String cimport *
from IsotopeDistribution cimport *
from Types cimport *
from Peak1D cimport *
from MassTrace cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>" namespace "OpenMS":

    cdef cppclass FLASHHelperClasses:
        # wrap-doc:
        #  Wrapper struct for all the structs needed by FLASHDeconv.
        #  Contains: PrecalculatedAveragine, MassFeature, IsobaricQuantities, LogMzPeak

        FLASHHelperClasses() except + nogil
        FLASHHelperClasses(FLASHHelperClasses &) except + nogil


cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>" namespace "OpenMS::FLASHHelperClasses":

    double getLogMz(double mz, bool positive) except + nogil  # wrap-attach:FLASHHelperClasses
    # wrap-doc:Calculate log mz from mz

    float getChargeMass(bool positive_ionization_mode) except + nogil  # wrap-attach:FLASHHelperClasses
    # wrap-doc:Get charge carrier mass (PROTON_MASS_U for positive, -PROTON_MASS_U for negative)


cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>" namespace "OpenMS":

    cdef cppclass MassFeature_FDHS "OpenMS::FLASHHelperClasses::MassFeature":
        # wrap-doc:
        #  Mass feature (Deconvolved masses in spectra are traced to generate mass features).
        #  Similar to LC-MS features but for deconvolved masses.

        # Constructors
        MassFeature_FDHS() except + nogil
        MassFeature_FDHS(MassFeature_FDHS &) except + nogil

        # Comparison operators
        bool operator<(MassFeature_FDHS& a) except + nogil
        bool operator>(MassFeature_FDHS& a) except + nogil
        bool operator==(MassFeature_FDHS& other) except + nogil

        # Public member variables
        UInt index  # wrap-doc:Feature index
        Kernel_MassTrace mt  # wrap-doc:The mass trace
        libcpp_vector[float] per_charge_intensity  # wrap-doc:Per charge intensities
        libcpp_vector[float] per_isotope_intensity  # wrap-doc:Per isotope intensities
        int iso_offset  # wrap-doc:Isotope offset between deconvolved masses and feature
        int scan_number  # wrap-doc:Representative scan number
        int min_scan_number  # wrap-doc:Minimum scan number
        int max_scan_number  # wrap-doc:Maximum scan number
        int rep_charge  # wrap-doc:Representative charge
        double avg_mass  # wrap-doc:Average mass
        int min_charge  # wrap-doc:Minimum charge
        int max_charge  # wrap-doc:Maximum charge
        int charge_count  # wrap-doc:Number of charges
        double isotope_score  # wrap-doc:Isotope score
        double qscore  # wrap-doc:Quality score
        double rep_mz  # wrap-doc:Representative m/z
        bool is_decoy  # wrap-doc:True if decoy
        UInt ms_level  # wrap-doc:MS level


    cdef cppclass PrecalAveragine "OpenMS::FLASHHelperClasses::PrecalculatedAveragine":
        # wrap-doc:
        #  Averagine patterns pre-calculated for speed up.
        #  Used for fast isotope cosine calculation.

        # Constructors
        PrecalAveragine() except + nogil
        PrecalAveragine(double min_mass, double max_mass, double delta, CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine) except + nogil
        PrecalAveragine(double min_mass, double max_mass, double delta, CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine, double decoy_iso_distance) except + nogil
        PrecalAveragine(PrecalAveragine &) except + nogil

        # Methods
        IsotopeDistribution get(double mass) except + nogil  # wrap-doc:Get isotope distribution for given mass
        size_t getMaxIsotopeIndex() except + nogil  # wrap-doc:Get max isotope index
        void setMaxIsotopeIndex(int index) except + nogil  # wrap-doc:Set max isotope index
        Size getLeftCountFromApex(double mass) except + nogil  # wrap-doc:Get isotope count left of apex
        Size getRightCountFromApex(double mass) except + nogil  # wrap-doc:Get isotope count right of apex
        Size getApexIndex(double mass) except + nogil  # wrap-doc:Get apex isotope index
        Size getLastIndex(double mass) except + nogil  # wrap-doc:Get last isotope index
        double getAverageMassDelta(double mass) except + nogil  # wrap-doc:Get mass diff between avg and mono
        double getMostAbundantMassDelta(double mass) except + nogil  # wrap-doc:Get mass diff between most abundant and mono
        double getSNRMultiplicationFactor(double mass) except + nogil  # wrap-doc:Get SNR multiplication factor


    cdef cppclass IsobaricQuantities "OpenMS::FLASHHelperClasses::IsobaricQuantities":
        # wrap-doc:
        #  Isobaric quantities from isobaric quantification.

        # Constructors
        IsobaricQuantities() except + nogil
        IsobaricQuantities(IsobaricQuantities &) except + nogil

        # Method
        bool empty() except + nogil  # wrap-doc:Returns true if no quantities stored

        # Public member variables
        int scan  # wrap-doc:Scan number
        double rt  # wrap-doc:Retention time
        double precursor_mz  # wrap-doc:Precursor m/z
        double precursor_mass  # wrap-doc:Precursor mass
        libcpp_vector[double] quantities  # wrap-doc:Isobaric quantities
        libcpp_vector[double] merged_quantities  # wrap-doc:Merged isobaric quantities


    cdef cppclass LogMzPeak "OpenMS::FLASHHelperClasses::LogMzPeak":
        # wrap-doc:
        #  Log transformed peak from original peak.
        #  Contains information such as charge, isotope index, and uncharged mass.

        # Constructors
        LogMzPeak() except + nogil
        LogMzPeak(LogMzPeak &) except + nogil
        LogMzPeak(Peak1D & peak, bool positive) except + nogil  # wrap-doc:Constructor from Peak1D

        # Methods
        double getUnchargedMass() except + nogil  # wrap-doc:Get uncharged mass (0 if no charge set)

        # Comparison operators
        bool operator<(LogMzPeak& a) except + nogil
        bool operator>(LogMzPeak& a) except + nogil
        bool operator==(LogMzPeak& other) except + nogil

        # Public member variables
        double mz  # wrap-doc:Original peak m/z
        float intensity  # wrap-doc:Original peak intensity
        double logMz  # wrap-doc:Log transformed m/z
        double mass  # wrap-doc:Determined mass after deconvolution (decharged, not monoisotopic)
        int abs_charge  # wrap-doc:Absolute charge
        bool is_positive  # wrap-doc:True if positive ionization mode
        int isotopeIndex  # wrap-doc:Isotope index
