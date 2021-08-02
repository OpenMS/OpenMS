from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/Base64.h>" namespace "OpenMS":
    
    cdef cppclass Base64 "OpenMS::Base64":

        Base64() nogil except + # wrap-doc:Class to encode and decode Base64, it supports two precisions 32 bit (float) and 64 bit (double).
        Base64(Base64 &) nogil except + # compiler

        void encode(libcpp_vector[ double ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except +#wrap-ignore
        void decode(const String & in_, ByteOrder from_byte_order, libcpp_vector[ double ] &out, bool zlib_compression) nogil except +#wrap-ignore
        void encode(libcpp_vector[ float ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except +#wrap-ignore
        void decode(const String & in_, ByteOrder from_byte_order, libcpp_vector[ float ] &out, bool zlib_compression) nogil except +#wrap-ignore

        void encode64(libcpp_vector[ double ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except + #wrap-ignore
        void decode64(const String & in_, ByteOrder from_byte_order, libcpp_vector[ double ] &out, bool zlib_compression) nogil except + #wrap-ignore
        void encode32(libcpp_vector[ float ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except + #wrap-ignore
        void decode32(const String & in_, ByteOrder from_byte_order, libcpp_vector[ float ] &out, bool zlib_compression) nogil except + #wrap-ignore

        void encodeIntegers(libcpp_vector[ int ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except + # wrap-doc:Encodes a vector of integer point numbers to a Base64 string
        void decodeIntegers(const String & in_, ByteOrder from_byte_order, libcpp_vector[ int ] &out, bool zlib_compression) nogil except + # wrap-doc:Decodes a Base64 string to a vector of integer numbers

        void encodeStrings(libcpp_vector[ String ] & in_, String &out, bool zlib_compression) nogil except + # wrap-doc:Encodes a vector of strings to a Base64 string
        void decodeStrings(const String & in_, libcpp_vector[ String ] &out, bool zlib_compression) nogil except + # wrap-doc:Decodes a Base64 string to a vector of (null-terminated) strings

        # void decodeSingleString(const String & in, QByteArray & base64_uncompressed, bool zlib_compression);

cdef extern from "<OpenMS/FORMAT/Base64.h>" namespace "OpenMS::Base64":
    cdef enum ByteOrder "OpenMS::Base64::ByteOrder":
        #wrap-attach:
        #    Base64
        BYTEORDER_BIGENDIAN
        BYTEORDER_LITTLEENDIAN
