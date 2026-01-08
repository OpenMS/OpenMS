from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from MSSpectrum cimport *
from FLASHHelperClasses cimport *
from PeakGroup cimport *
from Precursor cimport *
from Types cimport *


cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>" namespace "OpenMS":

    cdef cppclass DeconvolvedSpectrum:

        # default constructor
        DeconvolvedSpectrum() except + nogil
        # copy constructor
        DeconvolvedSpectrum(DeconvolvedSpectrum &) except + nogil
        # constructor with scan number
        DeconvolvedSpectrum(int scan_number) except + nogil

        # comparison operators
        bool operator<(DeconvolvedSpectrum & a) except + nogil
        bool operator>(DeconvolvedSpectrum & a) except + nogil
        bool operator==(DeconvolvedSpectrum & a) except + nogil

        # getters
        MSSpectrum getOriginalSpectrum() except + nogil
        PeakGroup getPrecursorPeakGroup() except + nogil
        int getPrecursorCharge() except + nogil
        Precursor getPrecursor() except + nogil
        double getCV() except + nogil
        double getCurrentMaxMass(double max_mass) except + nogil
        double getCurrentMinMass(double min_mass) except + nogil
        int getCurrentMaxAbsCharge(int max_abs_charge) except + nogil
        int getScanNumber() except + nogil
        int getPrecursorScanNumber() except + nogil
        bool isDecoy() except + nogil

        # setters
        void setPrecursor(Precursor & precursor) except + nogil
        void setPrecursorScanNumber(int scan_number) except + nogil
        void setPrecursorPeakGroup(PeakGroup & pg) except + nogil
        void setOriginalSpectrum(MSSpectrum & spec) except + nogil
        void setPeakGroups(libcpp_vector[PeakGroup] & x) except + nogil

        # vector-like interface for PeakGroups
        # wrap-iter-begin:__iter__(PeakGroup)
        # wrap-iter-end:__iter__(PeakGroup)
        Size size() except + nogil
        bool empty() except + nogil
        void clear() except + nogil
        void reserve(Size n) except + nogil
        void push_back(PeakGroup & pg) except + nogil
        void emplace_back(PeakGroup & pg) except + nogil
        void pop_back() except + nogil
        PeakGroup operator[](Size i) except + nogil
        PeakGroup back() except + nogil

        # sorting methods
        void sort() except + nogil
        void sortByQscore() except + nogil

        # conversion
        MSSpectrum toSpectrum(int to_charge, double tol, bool retain_undeconvolved) except + nogil
