from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *
# from SparseVector cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>" namespace "OpenMS":
    
    cdef cppclass BinnedSpectrum(MetaInfoInterface):
        BinnedSpectrum() nogil except +
        BinnedSpectrum(BinnedSpectrum) nogil except +
        BinnedSpectrum(Real size, UInt spread, MSSpectrum[Peak1D] ps) nogil except +
        bool operator==(BinnedSpectrum & rhs) nogil except +
        bool operator!=(BinnedSpectrum & rhs) nogil except +
        bool operator==(MSSpectrum[Peak1D] & rhs) nogil except +
        bool operator!=(MSSpectrum[Peak1D] & rhs) nogil except +
        double getBinSize() nogil except +
        UInt getBinSpread() nogil except +
        UInt getBinNumber() nogil except +
        UInt getFilledBinNumber() nogil except +
        ## SparseVector[ Real ]  getBins() nogil except +
        ## SparseVector[ Real ]  getBins() nogil except +
        ## const_bin_iterator begin() nogil except +
        ## const_bin_iterator end() nogil except +
        ## bin_iterator begin() nogil except +
        ## bin_iterator end() nogil except +
        void setBinSize(double s) nogil except +
        void setBinSpread(UInt s) nogil except +
        void setBinning() nogil except +
        bool checkCompliance(BinnedSpectrum & bs) nogil except +

        double getRT() nogil except +
        void   setRT(double) nogil except +
        unsigned int getMSLevel() nogil except +
        void setMSLevel(unsigned int) nogil except +

        libcpp_string getName() nogil except +
