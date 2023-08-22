from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from Types cimport *
from ParamValue cimport *
from String cimport *
from ParamEntry cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":
    
    cdef cppclass ParamNode "OpenMS::Param::ParamNode":
        ParamNode() except + nogil  # TODO
        ParamNode(ParamNode &) except + nogil 

        String name
        String description
        libcpp_vector[ ParamEntry ] entries
        libcpp_vector[ ParamNode ] nodes
        ParamNode(const String & n, const String & d) except + nogil 
        bool operator==(ParamNode & rhs) except + nogil 
        # EntryIterator findEntry(const String & name) except + nogil 
        # NodeIterator findNode(const String & name) except + nogil 
        ParamNode * findParentOf(const String & name) except + nogil 
        ParamEntry * findEntryRecursive(const String & name) except + nogil 
        void insert(ParamNode & node, const String & prefix) except + nogil 
        void insert(ParamEntry & entry, const String & prefix) except + nogil 
        Size size() except + nogil 
        String suffix(const String & key) except + nogil 

