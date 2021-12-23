from Types cimport *
from libcpp cimport bool
from IMSElement cimport *
from String cimport *
# from IMSAlphabetParser cimport *

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>" namespace "OpenMS::ims::IMSAlphabet":

    ctypedef IMSElement element_type
    ctypedef libcpp_vector[element_type] container
    ctypedef libcpp_vector[mass_type] masses_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>" namespace "OpenMS::ims":
    
    cdef cppclass IMSAlphabet "OpenMS::ims::IMSAlphabet":
        # wrap-doc:
                #   Holds an indexed list of bio-chemical elements.
                #   -----
                #   Presents an indexed list of bio-chemical elements of type (or derived from
                #   type) 'Element'. Due to indexed structure 'Alphabet' can be used similar
                #   to std::vector, for example to add a new element to 'Alphabet' function
                #   push_back(element_type) can be used. Elements or their properties (such
                #   as element's mass) can be accessed by index in a constant time. On the other
                #   hand accessing elements by their names takes linear time. Due to this and
                #   also the fact that 'Alphabet' is 'heavy-weighted' (consisting of
                #   'Element' -s or their derivatives where the depth of derivation as well is
                #   undefined resulting in possibly 'heavy' access operations) it is recommended
                #   not use 'Alphabet' directly in operations where fast access to
                #   'Element' 's properties is required. Instead consider to use
                #   'light-weighted' equivalents, such as 'Weights'
                #   -----
                #   :param map: MSExperiment to receive the identifications
                #   :param fmap: FeatureMap with PeptideIdentifications for the MSExperiment
                #   :param clear_ids: Reset peptide and protein identifications of each scan before annotating
                #   :param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z) 
                
        IMSAlphabet() nogil except +
        IMSAlphabet(IMSAlphabet &) nogil except +

        element_type  getElement(name_type & name) nogil except + # wrap-doc:Gets the element with 'index' and returns element with the given index in alphabet
        name_type  getName(size_type index) nogil except + # wrap-doc:Gets the symbol of the element with an 'index' in alphabet
        mass_type getMass(name_type & name) nogil except + # wrap-doc:Gets mono isotopic mass of the element with the symbol 'name' 
        mass_type getMass(size_type index) nogil except + # wrap-doc:Gets mass of the element with an 'index' in alphabet
        masses_type getMasses(size_type isotope_index) nogil except + # wrap-doc:Gets masses of elements isotopes given by 'isotope_index'
        masses_type getAverageMasses() nogil except + # wrap-doc:Gets average masses of elements
        bool hasName(name_type & name) nogil except + # wrap-doc:Returns true if there is an element with symbol 'name' in the alphabet, false - otherwise
        void push_back(name_type & name, mass_type value) nogil except + # wrap-doc:Adds a new element with 'name' and mass 'value'
        void push_back(element_type & element) nogil except + # wrap-doc:Adds a new 'element' to the alphabet
        void clear() nogil except + # wrap-doc:Clears the alphabet data
        void sortByNames() nogil except + # wrap-doc:Sorts the alphabet by names
        void sortByValues() nogil except + # wrap-doc:Sorts the alphabet by mass values
        void load(String & fname) nogil except + # wrap-doc:Loads the alphabet data from the file 'fname' using the default parser. If there is no file 'fname', throws an 'IOException'
        # POINTER # void load(libcpp_string & fname, IMSAlphabetParser[] * parser) nogil except +
        IMSAlphabet(libcpp_vector[IMSElement] & elements) nogil except +
        size_type size() nogil except +
        element_type  getElement(size_type index) nogil except + # wrap-doc:Gets the element with 'index'
        void setElement(name_type & name, mass_type mass, bool forced) nogil except + # wrap-doc:Overwrites an element in the alphabet with the 'name' with a new element constructed from the given 'name' and 'mass'
        bool erase(name_type & name) nogil except + # wrap-doc:Removes the element with 'name' from the alphabet

