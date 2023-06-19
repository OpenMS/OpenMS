from Types cimport *
from libcpp cimport bool
from IMSIsotopeDistribution cimport *
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>" namespace "OpenMS::ims::IMSElement":

    ctypedef libcpp_string name_type
    ctypedef IMSIsotopeDistribution isotopes_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>" namespace "OpenMS::ims":
    
    cdef cppclass IMSElement "OpenMS::ims::IMSElement":
        IMSElement() nogil except + # wrap-doc:Represents a chemical atom with name and isotope distribution
        IMSElement(IMSElement &) nogil except +
    
        # mass_type ELECTRON_MASS_IN_U read-only
        IMSElement(name_type & name, isotopes_type & isotopes) nogil except +
        IMSElement(name_type & name, mass_type mass) nogil except +
        IMSElement(name_type & name, nominal_mass_type nominal_mass) nogil except +
        name_type  getName() nogil except + # wrap-doc:Gets element's name
        void setName(name_type & name) nogil except + # wrap-doc:Sets element's name
        name_type  getSequence() nogil except + # wrap-doc:Gets element's sequence
        void setSequence(name_type & sequence) nogil except + # wrap-doc:Sets element's sequence
        nominal_mass_type getNominalMass() nogil except + # wrap-doc:Gets element's nominal mass
        mass_type getMass(size_type index) nogil except + # wrap-doc:Gets mass of element's isotope 'index'
        mass_type getAverageMass() nogil except + # wrap-doc:Gets element's average mass
        mass_type getIonMass(int electrons_number) nogil except + # wrap-doc:Gets ion mass of element. By default ion lacks 1 electron, but this can be changed by setting other 'electrons_number'
        IMSIsotopeDistribution  getIsotopeDistribution() nogil except + # wrap-doc:Gets element's isotope distribution
        void setIsotopeDistribution(IMSIsotopeDistribution & isotopes) nogil except + # wrap-doc:Sets element's isotope distribution
        bool operator==(IMSElement & element) nogil except +
        bool operator!=(IMSElement & element) nogil except +

