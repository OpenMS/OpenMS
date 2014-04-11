from Types cimport *
from Base64 cimport *

cdef extern from "<OpenMS/FORMAT/MSNumpressCoder.h>" namespace "OpenMS":

    cdef cppclass MSNumpressCoder:

        MSNumpressCoder() nogil except +
        MSNumpressCoder(MSNumpressCoder) nogil except +

        void encodeNP(libcpp_vector[double] in_, String & result,
                bool zlib_compression, NumpressConfig config) nogil except +

        void decodeNP(String in_, libcpp_vector[double] & out,
                bool zlib_compression, NumpressConfig config) nogil except +

cdef extern from "<OpenMS/FORMAT/MSNumpressCoder.h>" namespace "OpenMS::MSNumpressCoder":

    cdef enum NumpressCompression:
      # wrap-attach:
      #     MSNumpressCoder
      NONE,
      LINEAR,
      PIC,
      SLOF

    cdef cppclass NumpressConfig:
      # wrap-attach:
      #     MSNumpressCoder

      NumpressConfig() nogil except +
      NumpressConfig(NumpressConfig) nogil except +

      double numpressFixedPoint
      double numpressErrorTolerance
      NumpressCompression np_compression
      bool estimate_fixed_point

