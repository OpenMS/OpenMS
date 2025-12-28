from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from FLASHHelperClasses cimport *


cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>" namespace "OpenMS":

    cdef cppclass PeakGroup:

        # default constructor
        PeakGroup() except + nogil
        # copy constructor
        PeakGroup(PeakGroup &) except + nogil
        # constructor with charge range
        PeakGroup(int min_abs_charge, int max_abs_charge, bool is_positive) except + nogil

        # comparison operators
        bool operator<(PeakGroup & a) except + nogil
        bool operator>(PeakGroup & a) except + nogil
        bool operator==(PeakGroup & a) except + nogil

        # getters
        int getScanNumber() except + nogil
        double getMonoMass() except + nogil
        float getIntensity() except + nogil
        float getChargeSNR(int abs_charge) except + nogil
        float getChargeIsotopeCosine(int abs_charge) except + nogil
        float getChargeIntensity(int abs_charge) except + nogil
        float getIsotopeCosine() except + nogil
        float getPeakOccupancy() except + nogil
        int getRepAbsCharge() except + nogil
        double getQscore() except + nogil
        double getQscore2D() except + nogil
        float getSNR() except + nogil
        float getChargeScore() except + nogil
        float getAvgPPMError() except + nogil
        float getAvgDaError() except + nogil
        bool isPositive() except + nogil
        bool isTargeted() except + nogil
        float getQvalue() except + nogil
        double getIsotopeDaDistance() except + nogil
        int getMinNegativeIsotopeIndex() except + nogil
        UInt getIndex() except + nogil
        UInt getFeatureIndex() except + nogil

        # setters
        void setScanNumber(int scan_number) except + nogil
        void setChargeIsotopeCosine(int abs_charge, float cos) except + nogil
        void setAbsChargeRange(int min_abs_charge, int max_abs_charge) except + nogil
        void setIsotopeCosine(float cos) except + nogil
        void setRepAbsCharge(int max_snr_abs_charge) except + nogil
        void setMonoisotopicMass(double mono_mass) except + nogil
        void setQscore(double qscore) except + nogil
        void setChargeScore(float charge_score) except + nogil
        void setAvgPPMError(float error) except + nogil
        void setSNR(float snr) except + nogil
        void setChargeSNR(int abs_charge, float c_snr) except + nogil
        void setTargeted() except + nogil
        void setQvalue(double q) except + nogil
        void setIsotopeDaDistance(double d) except + nogil
        void setIndex(UInt i) except + nogil
        void setQscore2D(double fqscore) except + nogil
        void setFeatureIndex(UInt findex) except + nogil

        # vector-like operations
        Size size() except + nogil
        bool empty() except + nogil
        void reserve(Size n) except + nogil
        void sort() except + nogil
        void push_back(LogMzPeak & pg) except + nogil
        LogMzPeak operator[](Size i) except + nogil
