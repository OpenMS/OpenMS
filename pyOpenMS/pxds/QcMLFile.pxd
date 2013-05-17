from Types cimport *
from libcpp cimport bool
from ProgressLogger cimport *
# from XMLHandler cimport *
from XMLFile cimport *
from String cimport *
from UniqueIdGenerator cimport *

cdef extern from "<OpenMS/FORMAT/QcMLFile.h>" namespace "OpenMS::QcMLFile":
    
    cdef cppclass QualityParameter "OpenMS::QcMLFile::QualityParameter":
        QualityParameter() nogil except +
        QualityParameter(QualityParameter) nogil except +
        String name
        String id
        String value
        String cvRef
        String cvAcc
        String unitRef
        String unitAcc
        String flag
        bool operator==(QualityParameter & rhs) nogil except +
        bool operator<(QualityParameter & rhs) nogil except +
        bool operator>(QualityParameter & rhs) nogil except +
        String toXMLString(UInt indentation_level) nogil except +

