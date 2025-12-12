from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool
from Types cimport *
from String cimport *
from MSSpectrum cimport *
from FLASHHelperClasses cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>" namespace "OpenMS":

    cdef cppclass PeakGroup:
        # wrap-doc:
        #  Class describing a deconvolved mass.
        #  A mass contains multiple (LogMz) peaks of different charges and isotope indices.
        #  PeakGroup is the set of such peaks representing a single monoisotopic mass.

        # Constructors
        PeakGroup() except + nogil
        PeakGroup(PeakGroup &) except + nogil  # copy constructor
        PeakGroup(int min_abs_charge, int max_abs_charge, bool is_positive) except + nogil

        # Comparison operators
        bool operator<(PeakGroup & a) except + nogil
        bool operator>(PeakGroup & a) except + nogil
        bool operator==(PeakGroup & a) except + nogil

        # Core getters
        double getMonoMass() except + nogil  # wrap-doc:Returns the monoisotopic mass
        float getIntensity() except + nogil  # wrap-doc:Returns the summed intensity
        double getQscore() except + nogil  # wrap-doc:Returns the quality score (0-1)
        double getQscore2D() except + nogil  # wrap-doc:Returns the 2D quality score incorporating feature-level information
        float getIsotopeCosine() except + nogil  # wrap-doc:Returns the isotope cosine score
        float getChargeScore() except + nogil  # wrap-doc:Returns the charge fit score
        float getSNR() except + nogil  # wrap-doc:Returns the signal-to-noise ratio
        int getRepAbsCharge() except + nogil  # wrap-doc:Returns the representative charge
        # getAbsChargeRange returns std::tuple which is not supported - use getMinAbsCharge/getMaxAbsCharge instead
        libcpp_vector[float] getIsotopeIntensities() except + nogil  # wrap-doc:Returns per-isotope intensities
        int getScanNumber() except + nogil  # wrap-doc:Returns the scan number
        UInt getIndex() except + nogil  # wrap-doc:Returns the peak group index
        UInt getFeatureIndex() except + nogil  # wrap-doc:Returns the feature index
        float getQvalue() except + nogil  # wrap-doc:Returns the q-value for FDR
        TargetDecoyType getTargetDecoyType() except + nogil  # wrap-doc:Returns target/decoy type
        bool isPositive() except + nogil  # wrap-doc:Returns true if positive ionization mode
        bool isTargeted() except + nogil  # wrap-doc:Returns true if this peak group was targeted
        float getPeakOccupancy() except + nogil  # wrap-doc:Returns peak occupancy (0-1)
        float getAvgPPMError() except + nogil  # wrap-doc:Returns average ppm error
        float getAvgDaError() except + nogil  # wrap-doc:Returns average Da error
        # getRepMzRange and getMzRange return std::tuple which is not supported
        float getChargeIntensity(int abs_charge) except + nogil  # wrap-doc:Returns intensity for given charge
        float getChargeSNR(int abs_charge) except + nogil  # wrap-doc:Returns SNR for given charge
        float getChargeIsotopeCosine(int abs_charge) except + nogil  # wrap-doc:Returns isotope cosine for given charge
        double getIsotopeDaDistance() except + nogil  # wrap-doc:Returns distance between consecutive isotopes
        int getMinNegativeIsotopeIndex() except + nogil  # wrap-doc:Returns minimum negative isotope index
        libcpp_vector[float] getMassErrors(bool ppm) except + nogil  # wrap-doc:Returns mass errors per isotope

        # Setters
        void setScanNumber(int scan_number) except + nogil  # wrap-doc:Sets the scan number
        void setMonoisotopicMass(double mono_mass) except + nogil  # wrap-doc:Sets the monoisotopic mass
        void setIsotopeCosine(float cos) except + nogil  # wrap-doc:Sets the isotope cosine score
        void setChargeScore(float charge_score) except + nogil  # wrap-doc:Sets the charge fit score
        void setSNR(float snr) except + nogil  # wrap-doc:Sets the SNR
        void setQscore(double qscore) except + nogil  # wrap-doc:Sets the quality score
        void setQscore2D(double fqscore) except + nogil  # wrap-doc:Sets the 2D quality score
        void setQvalue(double q) except + nogil  # wrap-doc:Sets the q-value
        void setRepAbsCharge(int max_snr_abs_charge) except + nogil  # wrap-doc:Sets the representative charge
        # setAbsChargeRange is internal - charge range is set via constructor
        void setIndex(UInt i) except + nogil  # wrap-doc:Sets the peak group index
        void setFeatureIndex(UInt findex) except + nogil  # wrap-doc:Sets the feature index
        void setTargetDecoyType(TargetDecoyType index) except + nogil  # wrap-doc:Sets the target/decoy type
        void setTargeted() except + nogil  # wrap-doc:Marks this peak group as targeted
        void setIsotopeDaDistance(double d) except + nogil  # wrap-doc:Sets isotope distance
        void setChargeIsotopeCosine(int abs_charge, float cos) except + nogil  # wrap-doc:Sets isotope cosine for given charge
        void setChargeSNR(int abs_charge, float c_snr) except + nogil  # wrap-doc:Sets SNR for given charge
        void setAvgPPMError(float error) except + nogil  # wrap-doc:Sets average ppm error

        # Container operations
        Size size() except + nogil  # wrap-doc:Returns number of LogMzPeaks
        bool empty() except + nogil  # wrap-doc:Returns true if no peaks
        void push_back(LogMzPeak & pg) except + nogil  # wrap-doc:Adds a LogMzPeak
        void reserve(Size n) except + nogil  # wrap-doc:Reserves space for n peaks
        void sort() except + nogil  # wrap-doc:Sorts peaks by log m/z

        # Element access (read-only, const reference)
        LogMzPeak operator[](Size i) except + nogil  # wrap-upper-limit:size()

        # Iterators
        libcpp_vector[LogMzPeak].iterator begin() except + nogil  # wrap-iter-begin:__iter__(LogMzPeak)
        libcpp_vector[LogMzPeak].iterator end() except + nogil  # wrap-iter-end:__iter__(LogMzPeak)


cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>" namespace "OpenMS::PeakGroup":

    cdef enum TargetDecoyType:
        # wrap-attach:
        #  PeakGroup
        target,
        noise_decoy,
        signal_decoy
