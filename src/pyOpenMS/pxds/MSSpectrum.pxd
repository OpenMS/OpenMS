from libcpp.vector cimport vector as libcpp_vector
from Peak1D cimport *
from String cimport *
from RangeManager cimport *
from DataArrays cimport *
from SpectrumSettings cimport *

# this class has addons, see the ./addons folder (../addons/MSSpectrum.pyx)

cdef extern from "<OpenMS/KERNEL/MSSpectrum.h>" namespace "OpenMS":

    cdef cppclass MSSpectrum(SpectrumSettings, RangeManagerMzInt):
        # wrap-inherits:
        #  SpectrumSettings
        #  RangeManagerMzInt
        #
        # wrap-doc:
        #   The representation of a 1D spectrum.
        #   Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
        #   Iterations yields access to underlying peak objects but is slower
        #   Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
        #   See help(SpectrumSettings) for information about meta-information
        #   -----
        #   Usage:
        #     ms_level = spectrum.getMSLevel()
        #     rt = spectrum.getRT()
        #     mz, intensities = spectrum.get_peaks()
        #
        #   Usage:
        #     from pyopenms import *
        #
        #     spectrum = MSSpectrum()
        #     spectrum.setDriftTime(25) # 25 ms
        #     spectrum.setRT(205.2) # 205.2 s
        #     spectrum.setMSLevel(3) # MS3
        #     p = Precursor()
        #     p.setIsolationWindowLowerOffset(1.5)
        #     p.setIsolationWindowUpperOffset(1.5)
        #     p.setMZ(600) # isolation at 600 +/- 1.5 Th
        #     p.setActivationEnergy(40) # 40 eV
        #     p.setCharge(4) # 4+ ion
        #     spectrum.setPrecursors( [p] )
        #
        #     # Add raw data to spectrum
        #     spectrum.set_peaks( ([401.5], [900]) )
        #
        #     # Additional data arrays / peak annotations
        #     fda = FloatDataArray()
        #     fda.setName("Signal to Noise Array")
        #     fda.push_back(15)
        #     sda = StringDataArray()
        #     sda.setName("Peak annotation")
        #     sda.push_back("y15++")
        #     spectrum.setFloatDataArrays( [fda] )
        #     spectrum.setStringDataArrays( [sda] )
        #
        #     # Add spectrum to MSExperiment
        #     exp = MSExperiment()
        #     exp.addSpectrum(spectrum)
        #
        #     # Add second spectrum and store as mzML file
        #     spectrum2 = MSSpectrum()
        #     spectrum2.set_peaks( ([1, 2], [1, 2]) )
        #     exp.addSpectrum(spectrum2)
        #
        #     MzMLFile().store("testfile.mzML", exp)
        #   -----

        MSSpectrum() nogil except +
        MSSpectrum(MSSpectrum &) nogil except +
        double getRT() nogil except + # wrap-doc:Returns the absolute retention time (in seconds)
        void setRT(double) nogil except +  # wrap-doc:Sets the absolute retention time (in seconds)
        double getDriftTime() nogil except + # wrap-doc:Returns the drift time (-1 if not set)
        void setDriftTime(double) nogil except + # wrap-doc:Sets the drift time (-1 if not set)
        unsigned int getMSLevel() nogil except + # wrap-doc:Returns the MS level
        void setMSLevel(unsigned int) nogil except + # wrap-doc:Sets the MS level

        String getName() nogil except +
        void setName(String) nogil except +

        Size size() nogil except + # wrap-doc:Returns the number of peaks in the spectrum
        void reserve(size_t n) nogil except + 

        Peak1D& operator[](int) nogil except + # wrap-upper-limit:size()

        void updateRanges() nogil except +
        void clear(bool clear_meta_data) nogil except + # wrap-doc:Clears all data (and meta data if clear_meta_data is true)
        void push_back(Peak1D)  nogil except + # wrap-doc:Append a peak

        bool isSorted() nogil except + # wrap-doc:Returns true if the spectrum is sorte by m/z

        int findNearest(double mz) nogil except + # wrap-doc:Returns the index of the closest peak in m/z
        int findNearest(double mz, double tolerance) nogil except + # wrap-doc:Returns the index of the closest peak in the provided +/- m/z tolerance window (-1 if none match)
        int findNearest(double mz, double tolerance_left, double tolerance_right) nogil except + # wrap-doc:Returns the index of the closest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)
        int findHighestInWindow(double mz, double tolerance_left, double tolerance_right) nogil except + # wrap-doc:Returns the index of the highest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)

        MSSpectrum select(libcpp_vector[ size_t ] & indices) nogil except + # wrap-doc:Subset the spectrum by indices. Also applies to associated data arrays if present.

        void assign(libcpp_vector[Peak1D].iterator, libcpp_vector[Peak1D].iterator) nogil except + # wrap-ignore
        libcpp_vector[Peak1D].iterator begin() nogil except +  # wrap-iter-begin:__iter__(Peak1D)
        libcpp_vector[Peak1D].iterator end()   nogil except +  # wrap-iter-end:__iter__(Peak1D)

        double calculateTIC() nogil except + # wrap-doc:Returns the total ion current (=sum) of peak intensities in the spectrum

        bool operator==(MSSpectrum) nogil except +
        bool operator!=(MSSpectrum) nogil except +

        void sortByIntensity(bool reverse) nogil except +
        void sortByPosition() nogil except +

        libcpp_vector[FloatDataArray] getFloatDataArrays() nogil except + # wrap-doc:Returns the additional float data arrays to store e.g. meta data
        libcpp_vector[IntegerDataArray] getIntegerDataArrays() nogil except + # wrap-doc:Returns the additional int data arrays to store e.g. meta data
        libcpp_vector[StringDataArray] getStringDataArrays() nogil except + # wrap-doc:Returns the additional string data arrays to store e.g. meta data

        void setFloatDataArrays(libcpp_vector[FloatDataArray] fda) nogil except + # wrap-doc:Sets the additional float data arrays to store e.g. meta data
        void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida) nogil except + # wrap-doc:Sets the additional int data arrays to store e.g. meta data
        void setStringDataArrays(libcpp_vector[StringDataArray] sda) nogil except + # wrap-doc:Sets the additional string data arrays to store e.g. meta data

