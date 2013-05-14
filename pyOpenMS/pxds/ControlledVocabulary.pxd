from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool

from String cimport *
from Map cimport *
from StringList cimport *
from CVMappings cimport *
from DataValue cimport *

cdef extern from "<OpenMS/FORMAT/ControlledVocabulary.h>" namespace "OpenMS":

    cdef cppclass ControlledVocabulary:

        ControlledVocabulary() nogil except +
        ControlledVocabulary(ControlledVocabulary & voc) nogil except +

        # Returns the CV name (set in the load method)
        String name() nogil except +

        void loadFromOBO(String name, String filename) nogil except +

        # Returns true if the term is in the CV. Returns false otherwise.
        bool exists(String id) nogil except +

        # Returns true if a term with the given name is in the CV. Returns false otherwise.
        bool hasTermWithName(String name) nogil except +

        CVTerm_ControlledVocabulary getTerm(String id) nogil except +

        CVTerm_ControlledVocabulary getTermByName(String name, String desc) nogil except +

        # returns all the terms stored in the CV
        # TODO OpenMS Map type
        # Map[String, CVTerm_ControlledVocabulary] getTerms() nogil except +

        void getAllChildTerms(libcpp_set[String] terms, String parent) nogil except +

        bool isChildOf(String child, String parent) nogil except +

cdef extern from "<OpenMS/FORMAT/ControlledVocabulary.h>" namespace "OpenMS::ControlledVocabulary":

    cdef cppclass CVTerm_ControlledVocabulary "OpenMS::ControlledVocabulary::CVTerm":

      String name  #< Text name
      String id  #< Identifier
      libcpp_set[String] parents  #< The parent IDs
      libcpp_set[String] children  #< The child IDs
      bool obsolete  #< Flag that indicates of the term is obsolete
      String description  #< Term description
      StringList synonyms  #< List of synonyms
      StringList unparsed  #< Unparsed lines from the definition file
      XRefType_CVTerm_ControlledVocabulary xref_type  #< xref value-type for the CV-term
      StringList xref_binary  #< xref binary-data-type for the CV-term (list of all allowed data value types for the current binary data array)
      libcpp_set[String] units  #< unit accession ids, defined by relationship has units

      #Default constructor
      CVTerm_ControlledVocabulary() nogil except +

      CVTerm_ControlledVocabulary(CVTerm_ControlledVocabulary rhs) nogil except +

      String toXMLString(String ref, String value) nogil except +
      String toXMLString(String ref, DataValue value) nogil except +
      String getXRefTypeName(XRefType_CVTerm_ControlledVocabulary type) nogil except +

cdef extern from "<OpenMS/FORMAT/ControlledVocabulary.h>" namespace "OpenMS::ControlledVocabulary::CVTerm":

    # define xsd types allowed in cv term to specify their value-type
    cdef enum XRefType_CVTerm_ControlledVocabulary "OpenMS::ControlledVocabulary::CVTerm::XRefType":
        XSD_STRING = 0, # xsd:string A string
        XSD_INTEGER, # xsd:integer Any integer
        XSD_DECIMAL, # xsd:decimal Any real number
        XSD_NEGATIVE_INTEGER, # xsd:negativeInteger Any negative integer
        XSD_POSITIVE_INTEGER, # xsd:positiveInteger Any integer ] 0
        XSD_NON_NEGATIVE_INTEGER, # xsd:nonNegativeInteger Any integer ]= 0
        XSD_NON_POSITIVE_INTEGER, # xsd:nonPositiveInteger Any integer < 0
        XSD_BOOLEAN, # xsd:boolean True or false
        XSD_DATE, # xsd:date An XML-Schema date
        XSD_ANYURI, # xsd:anyURI uniform resource identifier
        NONE

