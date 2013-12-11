from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":

  cdef cppclass OSBinaryDataArray:
        OSBinaryDataArray()
        OSBinaryDataArray(OSBinaryDataArray)
        libcpp_vector[double] data
      
  ctypedef shared_ptr[OSBinaryDataArray] OSBinaryDataArrayPtr

  # See addons/OSSpectrum.pyx
  cdef cppclass OSSpectrum:
        OSSpectrum()
        OSSpectrum(OSSpectrum)
        OSBinaryDataArrayPtr getMZArray() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray() #wrap-ignore
        void setMZArray(OSBinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(OSBinaryDataArrayPtr data) #wrap-ignore

  ctypedef shared_ptr[OSSpectrum] OSSpectrumPtr

  # See addons/OSChromatogram.pyx
  cdef cppclass OSChromatogram:
        OSChromatogram()
        OSChromatogram(OSChromatogram)
        OSBinaryDataArrayPtr getTimeArray() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray() #wrap-ignore
        void setTimeArray(OSBinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(OSBinaryDataArrayPtr data) #wrap-ignore

  ctypedef shared_ptr[OSChromatogram] OSChromatogramPtr

