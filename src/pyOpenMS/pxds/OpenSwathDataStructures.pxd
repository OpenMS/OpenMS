from Types cimport *
from smart_ptr cimport shared_ptr

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":

  cdef cppclass OSBinaryDataArray:
        OSBinaryDataArray() nogil except +
        OSBinaryDataArray(OSBinaryDataArray) nogil except +
        libcpp_vector[double] data
        libcpp_string description
      
  ctypedef shared_ptr[OSBinaryDataArray] OSBinaryDataArrayPtr

  # See ../addons/OSSpectrum.pyx
  cdef cppclass OSSpectrum:
        OSSpectrum() nogil except +
        OSSpectrum(OSSpectrum) nogil except +
        OSBinaryDataArrayPtr getMZArray() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray() #wrap-ignore
        # libcpp_vector[ BinaryDataArrayPtr ]  getDataArrays() nogil except +
        void setMZArray(OSBinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(OSBinaryDataArrayPtr data) #wrap-ignore

  # boost::shared_ptr<OpenSwath::Spectrum> OpenSwath::SpectrumPtr;
  ctypedef shared_ptr[OSSpectrum] OSSpectrumPtr

  # See ../addons/OSChromatogram.pyx
  cdef cppclass OSChromatogram:
        OSChromatogram() nogil except +
        OSChromatogram(OSChromatogram) nogil except +
        OSBinaryDataArrayPtr getTimeArray() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray() #wrap-ignore
        # libcpp_vector[ BinaryDataArrayPtr ]  getDataArrays() nogil except +
        void setTimeArray(OSBinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(OSBinaryDataArrayPtr data) #wrap-ignore

  # boost::shared_ptr<OpenSwath::Chromatogram> OpenSwath::ChromatogramPtr;
  ctypedef shared_ptr[OSChromatogram] OSChromatogramPtr

