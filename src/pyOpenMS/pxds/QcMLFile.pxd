from Types cimport *
from ProgressLogger cimport *
from XMLHandler cimport *
from XMLFile cimport *
from String cimport *
from StringList cimport *
from Attachment cimport *

cdef extern from "<OpenMS/FORMAT/QcMLFile.h>" namespace "OpenMS":
    
    cdef cppclass QcMLFile(XMLHandler,XMLFile,ProgressLogger) :
        # wrap-inherits:
        #  XMLHandler
        #  XMLFile
        #  ProgressLogger
        QcMLFile() nogil except +
        QcMLFile(QcMLFile) nogil except + #wrap-ignore
        String map2csv(libcpp_map[ String, libcpp_map[ String, String ] ] & cvs_table, const String & separator) nogil except + # wrap-ignore
        String exportIDstats(const String & filename) nogil except +
        void addRunQualityParameter(String r, QualityParameter qp) nogil except +
        void addRunAttachment(String r, Attachment at) nogil except +
        void addSetQualityParameter(String r, QualityParameter qp) nogil except +
        void addSetAttachment(String r, Attachment at) nogil except +
        void removeAttachment(String r, libcpp_vector[ String ] & ids, String at) nogil except +
        void removeAttachment(String r, String at) nogil except +
        void removeAllAttachments(String at) nogil except +
        void removeQualityParameter(String r, libcpp_vector[ String ] & ids) nogil except +
        void merge(QcMLFile & addendum, String setname) nogil except +
        void collectSetParameter(String setname, String qp, libcpp_vector[ String ] & ret) nogil except +
        String exportAttachment(String filename, String qpname) nogil except +
        void getRunNames(libcpp_vector[ String ] & ids) nogil except +
        bool existsRun(String filename) nogil except +
        bool existsSet(String filename) nogil except +
        void existsRunQualityParameter(String filename, String qpname, libcpp_vector[ String ] & ids) nogil except +
        void existsSetQualityParameter(String filename, String qpname, libcpp_vector[ String ] & ids) nogil except +
        void store(const String & filename) nogil except +
        void load(const String & filename) nogil except +

        void registerRun(String id_, String name) nogil except +
        void registerSet(String id_, String name, libcpp_set[ String ] & names) nogil except +
        String exportQP(String filename, String qpname) nogil except +
        String exportQPs(String filename, StringList qpnames) nogil except +
        void getRunIDs(libcpp_vector[ String ] & ids) nogil except +


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

