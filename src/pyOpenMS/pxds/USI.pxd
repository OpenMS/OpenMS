from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/USI.h>" namespace "OpenMS":

    cdef cppclass USI "OpenMS::USI":
        # wrap-doc:
        #  Utility class for handling Universal Spectrum Identifiers (USI).
        #
        #  USI format (PSI-MS MS:1003063):
        #    mzspec:<collection>:<ms_run>:<index_type>:<index>[:interpretation]
        #
        #  The optional interpretation part uses ProForma proteoform-ion notation.

        USI() except + nogil
        USI(USI &) except + nogil

        USI(const String& collection,
            const String& ms_run,
            IndexType index_type,
            const String& index,
            const String& interpretation) except + nogil

        USI(const String& usi_string) except + nogil

        bool isValid() except + nogil  # wrap-doc:Return True if required fields are set

        @staticmethod
        bool isValidUSI(const String& usi_string) except + nogil  # wrap-doc:Validate a USI string format

        String getCollection() except + nogil  # wrap-doc:Get the dataset/library identifier
        void setCollection(const String& collection) except + nogil  # wrap-doc:Set the dataset/library identifier

        String getMSRun() except + nogil  # wrap-doc:Get the MS run file name
        void setMSRun(const String& ms_run) except + nogil  # wrap-doc:Set the MS run file name

        IndexType getIndexType() except + nogil  # wrap-doc:Get the index type (scan/index/nativeId)
        void setIndexType(IndexType index_type) except + nogil  # wrap-doc:Set the index type

        String getIndex() except + nogil  # wrap-doc:Get the spectrum index value
        void setIndex(const String& index) except + nogil  # wrap-doc:Set the spectrum index value

        String getInterpretation() except + nogil  # wrap-doc:Get the optional ProForma interpretation
        void setInterpretation(const String& interpretation) except + nogil  # wrap-doc:Set the optional ProForma interpretation

        bool hasInterpretation() except + nogil  # wrap-doc:Return True if interpretation is present

        String toString() except + nogil  # wrap-doc:Convert this USI to its string representation (empty if invalid)
        bool fromString(const String& usi_string) except + nogil  # wrap-doc:Parse a USI string into this object

        @staticmethod
        String extractBasename(const String& filepath) except + nogil  # wrap-doc:Extract basename from file path/URI for use as ms_run

        @staticmethod
        String indexTypeToString(IndexType index_type) except + nogil  # wrap-doc:Convert index type enum to string

        @staticmethod
        IndexType indexTypeFromString(const String& type_string) except + nogil  # wrap-doc:Parse index type from string

        @staticmethod
        String getCVAccession() except + nogil  # wrap-doc:Get PSI-MS CV accession for USI (MS:1003063)

        @staticmethod
        String getCVName() except + nogil  # wrap-doc:Get PSI-MS CV name for USI


cdef extern from "<OpenMS/METADATA/USI.h>" namespace "OpenMS::USI":

    cdef enum class IndexType "OpenMS::USI::IndexType":
        # wrap-attach:
        #    USI
        SCAN
        INDEX
        NATIVEID
