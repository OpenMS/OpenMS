from Types cimport *
from Types cimport *
from DateTime cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/XMLHandler.h>" namespace "OpenMS::Internal":
    
    # inherits from xercesc DefaultHandler -> no wrapping of it xercesc classes
    cdef cppclass XMLHandler:
        # private
        XMLHandler() except + nogil  #wrap-ignore
        #  copy constructor of 'XMLHandler' is implicitly deleted because base class 'xercesc::DefaultHandler' has an inaccessible copy constructor public xercesc::DefaultHandler
        XMLHandler(XMLHandler &) except + nogil  # wrap-ignore

        # NAMESPACE # void fatalError(xercesc::SAXParseException & exception) except + nogil 
        # NAMESPACE # void error(xercesc::SAXParseException & exception) except + nogil 
        # NAMESPACE # void warning(xercesc::SAXParseException & exception) except + nogil 
        XMLHandler(const String & filename, const String & version) except + nogil 
        void reset() except + nogil 
        # TODO cdash might parse out "fatalError" statements and interpret them
        # as compilation failure...
        # void fatalError(ActionMode mode, const String & msg, UInt line, UInt column) except + nogil 
        void error(ActionMode mode, const String & msg, UInt line, UInt column) except + nogil 
        void warning(ActionMode mode, const String & msg, UInt line, UInt column) except + nogil 
        # POINTER # void characters(XMLCh *chars, XMLSize_t length) except + nogil 
        # NAMESPACE # # POINTER # void startElement(XMLCh *uri, XMLCh *localname, XMLCh *qname, xercesc::Attributes & attrs) except + nogil 
        # POINTER # void endElement(XMLCh *uri, XMLCh *localname, XMLCh *qname) except + nogil 
        # NAMESPACE # void writeTo(std::ostream & ) except + nogil 

cdef extern from "<OpenMS/FORMAT/HANDLERS/XMLHandler.h>" namespace "OpenMS::Internal::XMLHandler":
    cdef enum ActionMode "OpenMS::Internal::XMLHandler::ActionMode":
        #wrap-attach:
        #   XMLHandler
        LOAD
        STORE

