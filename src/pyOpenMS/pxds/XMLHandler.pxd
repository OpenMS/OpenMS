from Types cimport *
from Types cimport *
from DateTime cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/XMLHandler.h>" namespace "OpenMS::Internal":
    
    # inherits from xercesc DefaultHandler -> no wrapping of it xercesc classes
    cdef cppclass XMLHandler:
        XMLHandler() nogil except + #wrap-ignore
        XMLHandler(XMLHandler) nogil except + #wrap-ignore

        # NAMESPACE # void fatalError(xercesc::SAXParseException & exception) nogil except +
        # NAMESPACE # void error(xercesc::SAXParseException & exception) nogil except +
        # NAMESPACE # void warning(xercesc::SAXParseException & exception) nogil except +
        XMLHandler(const String & filename, const String & version) nogil except +
        void reset() nogil except +
        # TODO cdash might parse out "fatalError" statements and interpret them
        # as compilation failure...
        # void fatalError(ActionMode mode, const String & msg, UInt line, UInt column) nogil except +
        void error(ActionMode mode, const String & msg, UInt line, UInt column) nogil except +
        void warning(ActionMode mode, const String & msg, UInt line, UInt column) nogil except +
        # POINTER # void characters(XMLCh *chars, XMLSize_t length) nogil except +
        # NAMESPACE # # POINTER # void startElement(XMLCh *uri, XMLCh *localname, XMLCh *qname, xercesc::Attributes & attrs) nogil except +
        # POINTER # void endElement(XMLCh *uri, XMLCh *localname, XMLCh *qname) nogil except +
        # NAMESPACE # void writeTo(std::ostream & ) nogil except +
        String errorString() nogil except +

cdef extern from "<OpenMS/FORMAT/HANDLERS/XMLHandler.h>" namespace "OpenMS::Internal::XMLHandler":
    cdef enum ActionMode "OpenMS::Internal::XMLHandler::ActionMode":
        #wrap-attach:
        #    XMLHandler
        LOAD
        STORE

