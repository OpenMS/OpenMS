from Types cimport *
from libcpp cimport bool
from Types cimport *
from DataValue cimport *
from String cimport *
from Map cimport *
from ParamEntry cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":
    
    cdef cppclass ParamNode "OpenMS::Param::ParamNode":
        ParamNode() nogil except +
        ParamNode(ParamNode) nogil except + #wrap-ignore
        String name
        String description
        libcpp_vector[ ParamEntry ] entries
        libcpp_vector[ ParamNode ] nodes
        ParamNode(String & n, String & d) nogil except +
        bool operator==(ParamNode & rhs) nogil except +
        # EntryIterator findEntry(String & name) nogil except +
        # NodeIterator findNode(String & name) nogil except +
        # POINTER # ParamNode * findParentOf(String & name) nogil except +
        # POINTER # ParamEntry * findEntryRecursive(String & name) nogil except +
        void insert(ParamNode & node, String & prefix) nogil except +
        void insert(ParamEntry & entry, String & prefix) nogil except +
        Size size() nogil except +
        String suffix(String & key) nogil except +

