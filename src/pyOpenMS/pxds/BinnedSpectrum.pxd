from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *
# from SparseVector cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>" namespace "OpenMS":

    cdef cppclass BinnedSpectrum:

        BinnedSpectrum() nogil except +
        BinnedSpectrum(BinnedSpectrum) nogil except +
        BinnedSpectrum(float size, UInt spread, MSSpectrum[Peak1D] ps) nogil except +
        bool operator==(BinnedSpectrum & rhs) nogil except +
        bool operator!=(BinnedSpectrum & rhs) nogil except +
        bool operator==(MSSpectrum[Peak1D] & rhs) nogil except +
        bool operator!=(MSSpectrum[Peak1D] & rhs) nogil except +
        double getBinSize() nogil except +
        UInt getBinSpread() nogil except +
        UInt getBinNumber() nogil except +
        UInt getFilledBinNumber() nogil except +
        ## SparseVector[ float ]  getBins() nogil except +
        ## SparseVector[ float ]  getBins() nogil except +
        ## const_bin_iterator begin() nogil except +
        ## const_bin_iterator end() nogil except +
        ## bin_iterator begin() nogil except +
        ## bin_iterator end() nogil except +
        void setBinSize(double s) nogil except +
        void setBinSpread(UInt s) nogil except +
        void setBinning() nogil except +
        bool checkCompliance(BinnedSpectrum & bs) nogil except +
        MSSpectrum[Peak1D] getRawSpectrum() nogil except +

