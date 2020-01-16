from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/SourceFile.h>" namespace "OpenMS":

    cdef cppclass SourceFile:
        SourceFile() nogil except +
        SourceFile(SourceFile) nogil except +
        String getNameOfFile() nogil except +
        void setNameOfFile(String) nogil except +

        String getPathToFile() nogil except +
        void setPathToFile(String) nogil except +

        float getFileSize() nogil except +
        void setFileSize(float) nogil except +

        String getFileType() nogil except +
        void setFileType(String) nogil except +

        String getChecksum() nogil except +
        void setChecksum(String, ChecksumType) nogil except +

        ChecksumType getChecksumType() nogil except +

        String getNativeIDType() nogil except +
        void setNativeIDType(String) nogil except +

        String getNativeIDTypeAccession() nogil except +
        void setNativeIDTypeAccession(const String & accesssion) nogil except +

cdef extern from "<OpenMS/METADATA/SourceFile.h>" namespace "OpenMS::SourceFile":
    cdef enum ChecksumType:
           UNKNOWN_CHECKSUM, SHA1, MD5, SIZE_OF_CHECKSUMTYPE
