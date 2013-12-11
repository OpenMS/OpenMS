from smart_ptr cimport shared_ptr
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/INTERFACES/DataStructures.h>" namespace "OpenMS::Interfaces":

  cdef cppclass BinaryDataArray:
        BinaryDataArray()
        BinaryDataArray(BinaryDataArray)
        libcpp_vector[double] data
      
  ctypedef shared_ptr[BinaryDataArray] BinaryDataArrayPtr

  # See addons/Spectrum.pyx
  cdef cppclass Spectrum:
        Spectrum()
        Spectrum(Spectrum)
        BinaryDataArrayPtr getMZArray() #wrap-ignore
        BinaryDataArrayPtr getIntensityArray() #wrap-ignore
        void setMZArray(BinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(BinaryDataArrayPtr data) #wrap-ignore

  ctypedef shared_ptr[Spectrum] SpectrumPtr

  # See addons/Chromatogram.pyx
  cdef cppclass Chromatogram:
        Chromatogram()
        Chromatogram(Chromatogram)
        BinaryDataArrayPtr getTimeArray() #wrap-ignore
        BinaryDataArrayPtr getIntensityArray() #wrap-ignore
        void setTimeArray(BinaryDataArrayPtr data) #wrap-ignore
        void setIntensityArray(BinaryDataArrayPtr data) #wrap-ignore

  ctypedef shared_ptr[Chromatogram] ChromatogramPtr

