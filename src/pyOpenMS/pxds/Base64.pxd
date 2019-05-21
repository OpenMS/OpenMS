from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/Base64.h>" namespace "OpenMS":
    
    cdef cppclass Base64 "OpenMS::Base64":
        Base64() nogil except +
        Base64(Base64) nogil except + #wrap-ignore

        void encode(libcpp_vector[ double ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except +#wrap-ignore
        void decode(const String & in_, ByteOrder from_byte_order, libcpp_vector[ double ] &out, bool zlib_compression) nogil except +#wrap-ignore
        void encode(libcpp_vector[ float ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except +#wrap-ignore
        void decode(const String & in_, ByteOrder from_byte_order, libcpp_vector[ float ] &out, bool zlib_compression) nogil except +#wrap-ignore

        void encode64(libcpp_vector[ double ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except + #wrap-ignore
        void decode64(const String & in_, ByteOrder from_byte_order, libcpp_vector[ double ] &out, bool zlib_compression) nogil except + #wrap-ignore
        void encode32(libcpp_vector[ float ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except + #wrap-ignore
        void decode32(const String & in_, ByteOrder from_byte_order, libcpp_vector[ float ] &out, bool zlib_compression) nogil except + #wrap-ignore

        void encodeIntegers(libcpp_vector[ int ] & in_, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except +
        void decodeIntegers(const String & in_, ByteOrder from_byte_order, libcpp_vector[ int ] &out, bool zlib_compression) nogil except +

        void encodeStrings(libcpp_vector[ String ] & in_, String &out, bool zlib_compression) nogil except +
        void decodeStrings(const String & in_, libcpp_vector[ String ] &out, bool zlib_compression) nogil except +

        # void decodeSingleString(const String & in, QByteArray & base64_uncompressed, bool zlib_compression);

cdef extern from "<OpenMS/FORMAT/Base64.h>" namespace "OpenMS::Base64":
    cdef enum ByteOrder "OpenMS::Base64::ByteOrder":
        #wrap-attach:
        #    Base64
        BYTEORDER_BIGENDIAN
        BYTEORDER_LITTLEENDIAN
