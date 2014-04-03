from Types cimport *
from libcpp cimport bool
from libcpp.string cimport string as libcpp_string
from IMSElement cimport *
# from IMSAlphabetParser cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>" namespace "OpenMS::ims::IMSAlphabet":

    ctypedef IMSElement element_type
    ctypedef libcpp_vector[element_type] container
    ctypedef libcpp_vector[mass_type] masses_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>" namespace "OpenMS::ims":
    
    cdef cppclass IMSAlphabet "OpenMS::ims::IMSAlphabet":
        IMSAlphabet() nogil except +
        IMSAlphabet(IMSAlphabet) nogil except +
        element_type  getElement(name_type & name) nogil except +
        name_type  getName(size_type index) nogil except +
        mass_type getMass(name_type & name) nogil except +
        mass_type getMass(size_type index) nogil except +
        masses_type getMasses(size_type isotope_index) nogil except +
        masses_type getAverageMasses() nogil except +
        bool hasName(name_type & name) nogil except +
        void push_back(name_type & name, mass_type value) nogil except +
        void push_back(element_type & element) nogil except +
        void clear() nogil except +
        void sortByNames() nogil except +
        void sortByValues() nogil except +
        void load(libcpp_string & fname) nogil except +
        # POINTER # void load(libcpp_string & fname, IMSAlphabetParser[] * parser) nogil except +
        IMSAlphabet(libcpp_vector[IMSElement] & elements) nogil except +
        size_type size() nogil except +
        element_type  getElement(size_type index) nogil except +
        void setElement(name_type & name, mass_type mass, bool forced) nogil except +
        bool erase(name_type & name) nogil except +

