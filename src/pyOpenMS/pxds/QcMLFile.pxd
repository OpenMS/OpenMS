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
        # wrap-doc:
        #   File adapter for QcML files used to load and store QcML files
        #   -----
        #   This Class is supposed to internally collect the data for the qcML File

        QcMLFile() nogil except +
        # copy constructor of 'QcMLFile' is implicitly deleted because base class 'Internal::XMLHandler' has a deleted copy constructor public Internal::XMLHandler,
        QcMLFile(QcMLFile &) nogil except + # wrap-ignore
        String map2csv(libcpp_map[ String, libcpp_map[ String, String ] ] & cvs_table, const String & separator) nogil except + # wrap-ignore
        String exportIDstats(const String & filename) nogil except +
        void addRunQualityParameter(String r, QualityParameter qp) nogil except + # wrap-doc:Adds a QualityParameter to run by the name r
        void addRunAttachment(String r, Attachment at) nogil except + # wrap-doc:Adds a attachment to run by the name r
        void addSetQualityParameter(String r, QualityParameter qp) nogil except + # wrap-doc:Adds a QualityParameter to set by the name r
        void addSetAttachment(String r, Attachment at) nogil except + # wrap-doc:Adds a attachment to set by the name r
        void removeAttachment(String r, libcpp_vector[ String ] & ids, String at) nogil except + # wrap-doc:Removes attachments referencing an id given in ids, from run/set r. All attachments if no attachment name is given with at
        void removeAttachment(String r, String at) nogil except + # wrap-doc:Removes attachment with cv accession at from run/set r
        void removeAllAttachments(String at) nogil except + # wrap-doc:Removes attachment with cv accession at from  all runs/sets
        void removeQualityParameter(String r, libcpp_vector[ String ] & ids) nogil except + # wrap-doc:Removes QualityParameter going by one of the ID attributes given in ids
        void merge(QcMLFile & addendum, String setname) nogil except + # wrap-doc:Merges the given QCFile into this one
        void collectSetParameter(String setname, String qp, libcpp_vector[ String ] & ret) nogil except + # wrap-doc:Collects the values of given QPs (as CVid) of the given set
        String exportAttachment(String filename, String qpname) nogil except + # wrap-doc:Returns a String of a tab separated rows if found empty string else from run/set by the name filename of the qualityparameter by the name qpname
        void getRunNames(libcpp_vector[ String ] & ids) nogil except + # wrap-doc:Gives the names of the registered runs in the vector ids
        bool existsRun(String filename) nogil except + # wrap-doc:Returns true if the given run id is present in this file, if checkname is true it also checks the names
        bool existsSet(String filename) nogil except + # wrap-doc:Returns true if the given set id is present in this file, if checkname is true it also checks the names
        void existsRunQualityParameter(String filename, String qpname, libcpp_vector[ String ] & ids) nogil except + # wrap-doc:Returns the ids of the parameter name given if found in given run empty else
        void existsSetQualityParameter(String filename, String qpname, libcpp_vector[ String ] & ids) nogil except + # wrap-doc:Returns the ids of the parameter name given if found in given set, empty else
        void store(const String & filename) nogil except + # wrap-doc:Store the qcML file
        void load(const String & filename) nogil except + # wrap-doc:Load a QCFile

        void registerRun(String id_, String name) nogil except + # wrap-doc:Registers a run in the qcml file with the respective mappings
        void registerSet(String id_, String name, libcpp_set[ String ] & names) nogil except + # wrap-doc:Registers a set in the qcml file with the respective mappings
        String exportQP(String filename, String qpname) nogil except + # wrap-doc:Returns a String value in quotation of a QualityParameter by the name qpname in run/set by the name filename
        String exportQPs(String filename, StringList qpnames) nogil except + # wrap-doc:Returns a String of a tab separated QualityParameter by the name qpname in run/set by the name filename
        void getRunIDs(libcpp_vector[ String ] & ids) nogil except + # wrap-doc:Gives the ids of the registered runs in the vector ids

cdef extern from "<OpenMS/FORMAT/QcMLFile.h>" namespace "OpenMS::QcMLFile":
    
    cdef cppclass QualityParameter "OpenMS::QcMLFile::QualityParameter":
        QualityParameter() nogil except +
        QualityParameter(QualityParameter &) nogil except + # compiler
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
