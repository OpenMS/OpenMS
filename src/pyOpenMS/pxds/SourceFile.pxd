from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/SourceFile.h>" namespace "OpenMS":

    cdef cppclass SourceFile:
        SourceFile() except + nogil  # wrap-doc:Description of a file location, used to store the origin of (meta) data
        SourceFile(SourceFile &) except + nogil 
        String getNameOfFile() except + nogil  # wrap-doc:Returns the file name
        void setNameOfFile(String) except + nogil  # wrap-doc:Sets the file name

        String getPathToFile() except + nogil  # wrap-doc:Returns the file path
        void setPathToFile(String) except + nogil  # wrap-doc:Sets the file path

        float getFileSize() except + nogil  # wrap-doc:Returns the file size in MB
        void setFileSize(float) except + nogil  # wrap-doc:Sets the file size in MB

        String getFileType() except + nogil  # wrap-doc:Returns the file type
        void setFileType(String) except + nogil  # wrap-doc:Sets the file type

        String getChecksum() except + nogil  # wrap-doc:Returns the file's checksum
        void setChecksum(String, ChecksumType) except + nogil  # wrap-doc:Sets the file's checksum

        ChecksumType getChecksumType() except + nogil  # wrap-doc:Returns the checksum type

        String getNativeIDType() except + nogil  # wrap-doc:Returns the native ID type of the spectra
        void setNativeIDType(String) except + nogil  # wrap-doc:Sets the native ID type of the spectra

        String getNativeIDTypeAccession() except + nogil  # wrap-doc:Returns the nativeID of the spectra
        void setNativeIDTypeAccession(const String & accesssion) except + nogil  # wrap-doc:Sets the native ID of the spectra

cdef extern from "<OpenMS/METADATA/SourceFile.h>" namespace "OpenMS::SourceFile":
    cdef enum ChecksumType:
           UNKNOWN_CHECKSUM, SHA1, MD5, SIZE_OF_CHECKSUMTYPE
