from Types cimport *
from String cimport *
from Base64 cimport *

cdef extern from "<OpenMS/FORMAT/MSNumpressCoder.h>" namespace "OpenMS":

    cdef cppclass MSNumpressCoder:

        MSNumpressCoder() nogil except +
        MSNumpressCoder(MSNumpressCoder) nogil except +

        void encodeNP(libcpp_vector[double] in_,
                      String & result,
                      bool zlib_compression,
                      NumpressConfig config) nogil except +

        void decodeNP(const String& in_,
                     libcpp_vector[double] & out,
                     bool zlib_compression,
                     NumpressConfig config) nogil except +

        void encodeNPRaw(libcpp_vector[ double ] in_,
                         String & result, 
                         NumpressConfig config) nogil except +

        void decodeNPRaw(const String& in_,
                         libcpp_vector[ double ] & out,
                         NumpressConfig config) nogil except +

cdef extern from "<OpenMS/FORMAT/MSNumpressCoder.h>" namespace "OpenMS::MSNumpressCoder":

    cdef enum NumpressCompression:
      # wrap-attach:
      #     MSNumpressCoder
      NONE,
      LINEAR,
      PIC,
      SLOF,
      SIZE_OF_NUMPRESSCOMPRESSION

    cdef cppclass NumpressConfig:

      NumpressConfig() nogil except +
      NumpressConfig(NumpressConfig) nogil except +

      double numpressFixedPoint # fixed point for numpress algorithms
      double numpressErrorTolerance # check error tolerance after encoding, guarantee abs(1.0-(encoded/decoded)) <= this, 0=do not guarantee anything
      NumpressCompression np_compression # which compression schema to use
      bool estimate_fixed_point # whether to estimate the fixed point or use the one proved with numpressFixedPoint
      double linear_fp_mass_acc # desired mass accuracy for linear encoding (-1 no effect, use 0.0001 for 0.2 ppm accuracy @ 500 m/z)

      void setCompression(const String & compression) nogil except +

