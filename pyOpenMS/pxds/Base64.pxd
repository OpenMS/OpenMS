from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/Base64.h>" namespace "OpenMS":
    
    cdef cppclass Base64 "OpenMS::Base64":
        Base64() nogil except +
        Base64(Base64) nogil except + #wrap-ignore
        # TEMPLATE # void encode(libcpp_vector[ FromType ] &in, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except +
        # TEMPLATE # void decode(String &in, ByteOrder from_byte_order, libcpp_vector[ ToType ] &out, bool zlib_compression) nogil except +
        # TEMPLATE # void encodeIntegers(libcpp_vector[ FromType ] &in, ByteOrder to_byte_order, String &out, bool zlib_compression) nogil except +
        # TEMPLATE # void decodeIntegers(String &in, ByteOrder from_byte_order, libcpp_vector[ ToType ] &out, bool zlib_compression) nogil except +
        void encodeStrings(libcpp_vector[ String ] & in_, String &out, bool zlib_compression) nogil except +
        void decodeStrings(String & in_, libcpp_vector[ String ] &out, bool zlib_compression) nogil except +

cdef extern from "<OpenMS/FORMAT/Base64.h>" namespace "OpenMS::Base64":
    cdef enum ByteOrder "OpenMS::Base64::ByteOrder":
        #wrap-attach:
        #    Base64
        BYTEORDER_BIGENDIAN
        BYTEORDER_LITTLEENDIAN

