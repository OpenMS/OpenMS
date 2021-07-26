from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/INTERFACES/DataStructures.h>" namespace "OpenMS::Interfaces":

  cdef cppclass BinaryDataArray:
        # here we misuse wrap-instances for renaming the instance, wrap-as is not supported for
        # classes, only for methods
        #
        # wrap-instances:
        #   _Interfaces_BinaryDataArray := BinaryDataArray
        BinaryDataArray() nogil except +
        BinaryDataArray(BinaryDataArray &) nogil except +
        libcpp_vector[double] data

  ctypedef shared_ptr[BinaryDataArray] BinaryDataArrayPtr

  # See addons/Spectrum.pyx
  cdef cppclass Spectrum:
        # here we misuse wrap-instances for renaming the instance, wrap-as is not supported for
        # classes, only for methods
        #
        # wrap-instances:
        #   _Interfaces_Spectrum := Spectrum
        Spectrum() nogil except +
        Spectrum(Spectrum &) nogil except +
        BinaryDataArrayPtr getMZArray() nogil except + #wrap-ignore
        BinaryDataArrayPtr getIntensityArray() nogil except + #wrap-ignore
        void setMZArray(BinaryDataArrayPtr data) nogil except + #wrap-ignore
        void setIntensityArray(BinaryDataArrayPtr data) nogil except + #wrap-ignore

  ctypedef shared_ptr[Spectrum] SpectrumPtr

  # See addons/Chromatogram.pyx
  cdef cppclass Chromatogram:
        # here we misuse wrap-instances for renaming the instance, wrap-as is not supported for
        # classes, only for methods
        #
        # wrap-instances:
        #   _Interfaces_Chromatogram := Chromatogram
        Chromatogram() nogil except +
        Chromatogram(Chromatogram &) nogil except +
        BinaryDataArrayPtr getTimeArray() nogil except + #wrap-ignore
        BinaryDataArrayPtr getIntensityArray() nogil except + #wrap-ignore
        void setTimeArray(BinaryDataArrayPtr data) nogil except + #wrap-ignore
        void setIntensityArray(BinaryDataArrayPtr data) nogil except + #wrap-ignore

  ctypedef shared_ptr[Chromatogram] ChromatogramPtr

