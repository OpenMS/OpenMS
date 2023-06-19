from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/SourceFile.h>" namespace "OpenMS":

    cdef cppclass SourceFile:
        SourceFile() nogil except + # wrap-doc:Description of a file location, used to store the origin of (meta) data
        SourceFile(SourceFile &) nogil except +
        String getNameOfFile() nogil except + # wrap-doc:Returns the file name
        void setNameOfFile(String) nogil except + # wrap-doc:Sets the file name

        String getPathToFile() nogil except + # wrap-doc:Returns the file path
        void setPathToFile(String) nogil except + # wrap-doc:Sets the file path

        float getFileSize() nogil except + # wrap-doc:Returns the file size in MB
        void setFileSize(float) nogil except + # wrap-doc:Sets the file size in MB

        String getFileType() nogil except + # wrap-doc:Returns the file type
        void setFileType(String) nogil except + # wrap-doc:Sets the file type

        String getChecksum() nogil except + # wrap-doc:Returns the file's checksum
        void setChecksum(String, ChecksumType) nogil except + # wrap-doc:Sets the file's checksum

        ChecksumType getChecksumType() nogil except + # wrap-doc:Returns the checksum type

        String getNativeIDType() nogil except + # wrap-doc:Returns the native ID type of the spectra
        void setNativeIDType(String) nogil except + # wrap-doc:Sets the native ID type of the spectra

        String getNativeIDTypeAccession() nogil except + # wrap-doc:Returns the nativeID of the spectra
        void setNativeIDTypeAccession(const String & accesssion) nogil except + # wrap-doc:Sets the native ID of the spectra

cdef extern from "<OpenMS/METADATA/SourceFile.h>" namespace "OpenMS::SourceFile":
    cdef enum ChecksumType:
           UNKNOWN_CHECKSUM, SHA1, MD5, SIZE_OF_CHECKSUMTYPE
