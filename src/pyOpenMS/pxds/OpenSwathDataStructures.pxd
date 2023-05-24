from Types cimport *
from smart_ptr cimport shared_ptr

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>" namespace "OpenSwath":

  # See ../addons/OSBinaryDataArray.pyx
  cdef cppclass OSBinaryDataArray:
        OSBinaryDataArray() nogil except +
        OSBinaryDataArray(OSBinaryDataArray &) nogil except + # compiler
        libcpp_vector[double] data
        libcpp_string description

        libcpp_vector[double] getData(self) #wrap-ignore wrap-doc:Access to a copy of the underlying data using a numpy array
        libcpp_vector[double] getData_mv(self) #wrap-ignore wrap-doc:Access to the underlying data using a memory view
      
  ctypedef shared_ptr[OSBinaryDataArray] OSBinaryDataArrayPtr

  # See ../addons/OSSpectrum.pyx
  cdef cppclass OSSpectrum:
        OSSpectrum() nogil except +
        OSSpectrum(OSSpectrum &) nogil except + # compiler
        
        # Obtain a copy of the underlying data
        OSBinaryDataArrayPtr getMZArray() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray() #wrap-ignore
        OSBinaryDataArrayPtr getDriftTimeArray() #wrap-ignore

        # Obtain access to a memory view of the underlying data
        OSBinaryDataArrayPtr getMZArray_mv() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray_mv() #wrap-ignore
        OSBinaryDataArrayPtr getDriftTimeArray_mv() #wrap-ignore

        libcpp_vector[ OSBinaryDataArrayPtr ] getDataArrays() #wrap-ignore
        void setDataArrays( libcpp_vector[ OSBinaryDataArrayPtr ]) #wrap-ignore

        void setMZArray(OSBinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(OSBinaryDataArrayPtr data) #wrap-ignore

  # boost::shared_ptr<OpenSwath::Spectrum> OpenSwath::SpectrumPtr;
  ctypedef shared_ptr[OSSpectrum] OSSpectrumPtr

  # See ../addons/OSChromatogram.pyx
  cdef cppclass OSChromatogram:
        OSChromatogram() nogil except +
        OSChromatogram(OSChromatogram &) nogil except + # compiler

        # Obtain a copy of the underlying data
        OSBinaryDataArrayPtr getTimeArray() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray() #wrap-ignore

        # Obtain access to a memory view of the underlying data
        OSBinaryDataArrayPtr getTimeArray_mv() #wrap-ignore
        OSBinaryDataArrayPtr getIntensityArray_mv() #wrap-ignore

        libcpp_vector[ OSBinaryDataArrayPtr ] getDataArrays() #wrap-ignore
        void setDataArrays( libcpp_vector[ OSBinaryDataArrayPtr ]) #wrap-ignore

        void setTimeArray(OSBinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(OSBinaryDataArrayPtr data) #wrap-ignore

  # boost::shared_ptr<OpenSwath::Chromatogram> OpenSwath::ChromatogramPtr;
  ctypedef shared_ptr[OSChromatogram] OSChromatogramPtr

