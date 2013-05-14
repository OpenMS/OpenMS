from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/METADATA/SourceFile.h>" namespace "OpenMS":

    cdef cppclass SourceFile:
        SourceFile() nogil except +
        SourceFile(SourceFile) nogil except +
        libcpp_string getNameOfFile() nogil except +
        void setNameOfFile(libcpp_string) nogil except +

        libcpp_string getPathToFile() nogil except +
        void setPathToFile(libcpp_string) nogil except +

        float getFileSize() nogil except +
        void setFileSize(float) nogil except +

        libcpp_string getFileType() nogil except +
        void setFileType(libcpp_string) nogil except +

        libcpp_string getChecksum() nogil except +
        void setChecksum(libcpp_string, ChecksumType) nogil except +

        ChecksumType getChecksumType() nogil except +

        libcpp_string getNativeIDType() nogil except +
        void setNativeIDType(libcpp_string) nogil except +


cdef extern from "<OpenMS/METADATA/SourceFile.h>" namespace "OpenMS::SourceFile":
    cdef enum ChecksumType:
           UNKNOWN_CHECKSUM, SHA1, MD5, SIZE_OF_CHECKSUMTYPE
