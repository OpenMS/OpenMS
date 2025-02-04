from Types cimport *
from smart_ptr cimport shared_ptr

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":

  cdef cppclass OSBinaryDataArray:
        OSBinaryDataArray() except + nogil 
        OSBinaryDataArray(OSBinaryDataArray &) except + nogil  # compiler
        libcpp_vector[double] data
        libcpp_string description
      
  ctypedef shared_ptr[OSBinaryDataArray] OSBinaryDataArrayPtr

  # See ../addons/OSSpectrum.pyx
  cdef cppclass OSSpectrum:
        OSSpectrum() except + nogil 
        OSSpectrum(OSSpectrum &) except + nogil  # compiler
        OSBinaryDataArrayPtr getMZArray() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray() #wrap-ignore
        # libcpp_vector[ BinaryDataArrayPtr ]  getDataArrays() except + nogil 
        void setMZArray(OSBinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(OSBinaryDataArrayPtr data) #wrap-ignore

  # boost::shared_ptr<OpenSwath::Spectrum> OpenSwath::SpectrumPtr;
  ctypedef shared_ptr[OSSpectrum] OSSpectrumPtr

  # See ../addons/OSChromatogram.pyx
  cdef cppclass OSChromatogram:
        OSChromatogram() except + nogil 
        OSChromatogram(OSChromatogram &) except + nogil  # compiler
        OSBinaryDataArrayPtr getTimeArray() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray() #wrap-ignore
        # libcpp_vector[ BinaryDataArrayPtr ]  getDataArrays() except + nogil 
        void setTimeArray(OSBinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(OSBinaryDataArrayPtr data) #wrap-ignore

  # boost::shared_ptr<OpenSwath::Chromatogram> OpenSwath::ChromatogramPtr;
  ctypedef shared_ptr[OSChromatogram] OSChromatogramPtr

