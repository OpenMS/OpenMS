from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool

from String cimport *
from StringList cimport *
from CVMappings cimport *
from DataValue cimport *

cdef extern from "<OpenMS/FORMAT/ControlledVocabulary.h>" namespace "OpenMS":

    cdef cppclass ControlledVocabulary:

        ControlledVocabulary() except + nogil 
        ControlledVocabulary(ControlledVocabulary &) except + nogil  # compiler

        # Returns the CV name (set in the load method)
        String name() except + nogil  # wrap-doc:Returns the CV name (set in the load method)

        void loadFromOBO(String name, String filename) except + nogil  # wrap-doc:Loads the CV from an OBO file

        # Returns true if the term is in the CV. Returns false otherwise.
        bool exists(String id) except + nogil  # wrap-doc:Returns true if the term is in the CV. Returns false otherwise.

        # Returns true if a term with the given name is in the CV. Returns false otherwise.
        bool hasTermWithName(String name) except + nogil  # wrap-doc:Returns true if a term with the given name is in the CV. Returns false otherwise

        CVTerm_ControlledVocabulary getTerm(String id) except + nogil  # wrap-doc:Returns a term specified by ID

        CVTerm_ControlledVocabulary getTermByName(String name, String desc) except + nogil  # wrap-doc:Returns a term specified by name

        # returns all the terms stored in the CV
        # TODO OpenMS Map type
        # Map[String, CVTerm_ControlledVocabulary] getTerms() except + nogil 

        void getAllChildTerms(libcpp_set[String] terms, String parent) except + nogil  # wrap-doc:Writes all child terms recursively into terms

        bool isChildOf(String child, String parent) except + nogil  # wrap-doc:Returns True if `child` is a child of `parent`

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
      CVTerm_ControlledVocabulary() except + nogil 
      CVTerm_ControlledVocabulary(CVTerm_ControlledVocabulary rhs) except + nogil 

      String toXMLString(String ref, String value) except + nogil  # wrap-doc:Get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
      String toXMLString(String ref, DataValue value) except + nogil  # wrap-doc:Get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
      String getXRefTypeName(XRefType_CVTerm_ControlledVocabulary type) except + nogil 
      bool isHigherBetterScore(CVTerm_ControlledVocabulary term) except + nogil 

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

