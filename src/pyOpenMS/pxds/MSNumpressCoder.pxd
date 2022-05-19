from Types cimport *
from String cimport *
from Base64 cimport *

cdef extern from "<OpenMS/FORMAT/MSNumpressCoder.h>" namespace "OpenMS":

    cdef cppclass MSNumpressCoder:

        MSNumpressCoder() nogil except +
        MSNumpressCoder(MSNumpressCoder &) nogil except + # compiler

        void encodeNP(libcpp_vector[double] in_,
                      String & result,
                      bool zlib_compression,
                      NumpressConfig config) nogil except +
          # wrap-doc:
                #   Encodes a vector of floating point numbers into a Base64 string using numpress
                #   -----
                #   This code is obtained from the proteowizard implementation
                #   ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp (adapted by Hannes Roest)
                #   -----
                #   This function will first apply the numpress encoding to the data, then
                #   encode the result in base64 (with optional zlib compression before
                #   base64 encoding)
                #   -----
                #   :note In case of error, result string is empty
                #   -----
                #   :param in: The vector of floating point numbers to be encoded
                #   :param result: The resulting string
                #   :param zlib_compression: Whether to apply zlib compression after numpress compression
                #   :param config: The numpress configuration defining the compression strategy

        void decodeNP(const String& in_,
                     libcpp_vector[double] & out,
                     bool zlib_compression,
                     NumpressConfig config) nogil except +
          # wrap-doc:
                #   Decodes a Base64 string to a vector of floating point numbers using numpress
                #   -----
                #   This code is obtained from the proteowizard implementation
                #   ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp (adapted by Hannes Roest)
                #   -----
                #   This function will first decode the input base64 string (with optional
                #   zlib decompression after decoding) and then apply numpress decoding to
                #   the data
                #   -----
                #   :param in: The base64 encoded string
                #   :param out: The resulting vector of doubles
                #   :param zlib_compression: Whether to apply zlib de-compression before numpress de-compression
                #   :param config: The numpress configuration defining the compression strategy
                #   :raises:
                #     Exception: ConversionError if the string cannot be converted

        void encodeNPRaw(libcpp_vector[ double ] in_,
                         String & result, 
                         NumpressConfig config) nogil except +
          # wrap-doc:
                #   Encode the data vector "in" to a raw byte array
                #   -----
                #   :note In case of error, "result" is given back unmodified
                #   :note The result is not a string but a raw byte array and may contain zero bytes
                #   -----
                #   This performs the raw numpress encoding on a set of data and does no
                #   Base64 encoding on the result. Therefore the result string is likely
                #   *unsafe* to handle and is a raw byte array.
                #   -----
                #   Please use the safe versions above unless you need access to the raw
                #   byte arrays
                #   -----
                #   :param in: The vector of floating point numbers to be encoded
                #   :param result: The resulting string
                #   :param config: The numpress configuration defining the compression strategy

        void decodeNPRaw(const String& in_,
                         libcpp_vector[ double ] & out,
                         NumpressConfig config) nogil except +
          # wrap-doc:
                #   Decode the raw byte array "in" to the result vector "out"
                #   -----
                #   :note The string in should *only* contain the data and _no_ extra
                #   null terminating byte
                #   -----
                #   This performs the raw numpress decoding on a raw byte array (not Base64
                #   encoded). Therefore the input string is likely *unsafe* to handle and is
                #   basically a byte container
                #   -----
                #   Please use the safe versions above unless you need access to the raw
                #   byte arrays
                #   -----
                #   :param in: The base64 encoded string
                #   :param out: The resulting vector of doubles
                #   :param config: The numpress configuration defining the compression strategy

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

      NumpressConfig() nogil except + # compiler
      NumpressConfig(NumpressConfig &) nogil except + # compiler

      double numpressFixedPoint # fixed point for numpress algorithms
      double numpressErrorTolerance # check error tolerance after encoding, guarantee abs(1.0-(encoded/decoded)) <= this, 0=do not guarantee anything
      NumpressCompression np_compression # which compression schema to use
      bool estimate_fixed_point # whether to estimate the fixed point or use the one proved with numpressFixedPoint
      double linear_fp_mass_acc # desired mass accuracy for linear encoding (-1 no effect, use 0.0001 for 0.2 ppm accuracy @ 500 m/z)

      void setCompression(const String & compression) nogil except +
