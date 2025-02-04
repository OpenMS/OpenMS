from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/INTERFACES/DataStructures.h>" namespace "OpenMS::Interfaces":

  cdef cppclass BinaryDataArray:
        # here we misuse wrap-instances for renaming the instance, wrap-as is not supported for
        # classes, only for methods
        #
        # wrap-instances:
        #  _Interfaces_BinaryDataArray := BinaryDataArray
        BinaryDataArray() except + nogil 
        BinaryDataArray(BinaryDataArray &) except + nogil 
        libcpp_vector[double] data

  ctypedef shared_ptr[BinaryDataArray] BinaryDataArrayPtr

  # See addons/Spectrum.pyx
  cdef cppclass Spectrum:
        # here we misuse wrap-instances for renaming the instance, wrap-as is not supported for
        # classes, only for methods
        #
        # wrap-instances:
        #  _Interfaces_Spectrum := Spectrum
        Spectrum() except + nogil 
        Spectrum(Spectrum &) except + nogil 
        BinaryDataArrayPtr getMZArray() except + nogil  #wrap-ignore
        BinaryDataArrayPtr getIntensityArray() except + nogil  #wrap-ignore
        void setMZArray(BinaryDataArrayPtr data) except + nogil  #wrap-ignore
        void setIntensityArray(BinaryDataArrayPtr data) except + nogil  #wrap-ignore

  ctypedef shared_ptr[Spectrum] SpectrumPtr

  # See addons/Chromatogram.pyx
  cdef cppclass Chromatogram:
        # here we misuse wrap-instances for renaming the instance, wrap-as is not supported for
        # classes, only for methods
        #
        # wrap-instances:
        #  _Interfaces_Chromatogram := Chromatogram
        Chromatogram() except + nogil 
        Chromatogram(Chromatogram &) except + nogil 
        BinaryDataArrayPtr getTimeArray() except + nogil  #wrap-ignore
        BinaryDataArrayPtr getIntensityArray() except + nogil  #wrap-ignore
        void setTimeArray(BinaryDataArrayPtr data) except + nogil  #wrap-ignore
        void setIntensityArray(BinaryDataArrayPtr data) except + nogil  #wrap-ignore

  ctypedef shared_ptr[Chromatogram] ChromatogramPtr

