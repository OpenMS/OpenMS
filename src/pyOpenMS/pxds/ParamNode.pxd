from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from Types cimport *
from ParamValue cimport *
from String cimport *
from ParamEntry cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":
    
    cdef cppclass ParamNode "OpenMS::Param::ParamNode":
        ParamNode() nogil except + # TODO
        ParamNode(ParamNode &) nogil except +

        String name
        String description
        libcpp_vector[ ParamEntry ] entries
        libcpp_vector[ ParamNode ] nodes
        ParamNode(const String & n, const String & d) nogil except +
        bool operator==(ParamNode & rhs) nogil except +
        # EntryIterator findEntry(const String & name) nogil except +
        # NodeIterator findNode(const String & name) nogil except +
        ParamNode * findParentOf(const String & name) nogil except +
        ParamEntry * findEntryRecursive(const String & name) nogil except +
        void insert(ParamNode & node, const String & prefix) nogil except +
        void insert(ParamEntry & entry, const String & prefix) nogil except +
        Size size() nogil except +
        String suffix(const String & key) nogil except +

