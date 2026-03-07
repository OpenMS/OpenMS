from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from MSSpectrum cimport *
from Precursor cimport *
from FLASHHelperClasses cimport *
from PeakGroup cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>" namespace "OpenMS":

    cdef cppclass DeconvolvedSpectrum:
        # wrap-doc:
        #  A class representing a deconvolved spectrum.
        #  DeconvolvedSpectrum consists of PeakGroup instances representing masses.
        #  For MSn n>1, a PeakGroup representing the precursor mass is also added.

        # Constructors
        DeconvolvedSpectrum() except + nogil
        DeconvolvedSpectrum(DeconvolvedSpectrum &) except + nogil  # copy constructor
        DeconvolvedSpectrum(int scan_number) except + nogil  # wrap-doc:Constructor with scan number

        # Comparison operators
        bool operator<(DeconvolvedSpectrum & a) except + nogil
        bool operator>(DeconvolvedSpectrum & a) except + nogil
        bool operator==(DeconvolvedSpectrum & a) except + nogil

        # Conversion
        MSSpectrum toSpectrum(int to_charge, double tol, bool retain_undeconvolved) except + nogil
        # wrap-doc:
        #  Convert DeconvolvedSpectrum to MSSpectrum.
        #  :param to_charge: The charge of each peak in output
        #  :param tol: The ppm tolerance
        #  :param retain_undeconvolved: If true, undeconvolved peaks are included

        # Getters
        MSSpectrum getOriginalSpectrum() except + nogil  # wrap-doc:Returns the original spectrum
        PeakGroup getPrecursorPeakGroup() except + nogil  # wrap-doc:Returns the precursor peak group (MSn, n>1)
        int getPrecursorCharge() except + nogil  # wrap-doc:Returns the precursor charge
        Precursor getPrecursor() except + nogil  # wrap-doc:Returns the precursor peak
        int getScanNumber() except + nogil  # wrap-doc:Returns the scan number
        int getPrecursorScanNumber() except + nogil  # wrap-doc:Returns the precursor scan number
        double getCurrentMaxMass(double max_mass) except + nogil  # wrap-doc:Returns the current max mass
        double getCurrentMinMass(double min_mass) except + nogil  # wrap-doc:Returns the current min mass
        int getCurrentMaxAbsCharge(int max_abs_charge) except + nogil  # wrap-doc:Returns the current max charge
        IsobaricQuantities getQuantities() except + nogil  # wrap-doc:Returns isobaric quantities
        bool isDecoy() except + nogil  # wrap-doc:Returns true if this is a decoy spectrum

        # Setters
        void setOriginalSpectrum(MSSpectrum & spec) except + nogil  # wrap-doc:Sets the original spectrum
        void setPrecursor(Precursor & precursor) except + nogil  # wrap-doc:Sets the precursor
        void setPrecursorScanNumber(int scan_number) except + nogil  # wrap-doc:Sets the precursor scan number
        void setPrecursorPeakGroup(PeakGroup & pg) except + nogil  # wrap-doc:Sets the precursor peak group
        void setPeakGroups(libcpp_vector[PeakGroup] & x) except + nogil  # wrap-doc:Sets peak groups
        void setQuantities(IsobaricQuantities & quantities) except + nogil  # wrap-doc:Sets isobaric quantities

        # Container operations
        Size size() except + nogil  # wrap-doc:Returns number of peak groups
        bool empty() except + nogil  # wrap-doc:Returns true if no peak groups
        void clear() except + nogil  # wrap-doc:Clears all peak groups
        void reserve(Size n) except + nogil  # wrap-doc:Reserves space for n peak groups
        void push_back(PeakGroup & pg) except + nogil  # wrap-doc:Adds a peak group
        void pop_back() except + nogil  # wrap-doc:Removes the last peak group

        # Element access
        PeakGroup& operator[](Size i) except + nogil  # wrap-upper-limit:size()

        # Iterators
        libcpp_vector[PeakGroup].iterator begin() except + nogil  # wrap-iter-begin:__iter__(PeakGroup)
        libcpp_vector[PeakGroup].iterator end() except + nogil  # wrap-iter-end:__iter__(PeakGroup)

        # Sorting
        void sort() except + nogil  # wrap-doc:Sorts peak groups by monoisotopic mass
        void sortByQscore() except + nogil  # wrap-doc:Sorts peak groups by Qscore
