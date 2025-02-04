from Types cimport *
from libcpp cimport bool
from IMSIsotopeDistribution cimport *
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>" namespace "OpenMS::ims::IMSElement":

    ctypedef libcpp_string name_type
    ctypedef IMSIsotopeDistribution isotopes_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>" namespace "OpenMS::ims":
    
    cdef cppclass IMSElement "OpenMS::ims::IMSElement":
        IMSElement() except + nogil  # wrap-doc:Represents a chemical atom with name and isotope distribution
        IMSElement(IMSElement &) except + nogil 
    
        # mass_type ELECTRON_MASS_IN_U read-only
        IMSElement(name_type & name, isotopes_type & isotopes) except + nogil 
        IMSElement(name_type & name, mass_type mass) except + nogil 
        IMSElement(name_type & name, nominal_mass_type nominal_mass) except + nogil 
        name_type  getName() except + nogil  # wrap-doc:Gets element's name
        void setName(name_type & name) except + nogil  # wrap-doc:Sets element's name
        name_type  getSequence() except + nogil  # wrap-doc:Gets element's sequence
        void setSequence(name_type & sequence) except + nogil  # wrap-doc:Sets element's sequence
        nominal_mass_type getNominalMass() except + nogil  # wrap-doc:Gets element's nominal mass
        mass_type getMass(size_type index) except + nogil  # wrap-doc:Gets mass of element's isotope 'index'
        mass_type getAverageMass() except + nogil  # wrap-doc:Gets element's average mass
        mass_type getIonMass(int electrons_number) except + nogil  # wrap-doc:Gets ion mass of element. By default ion lacks 1 electron, but this can be changed by setting other 'electrons_number'
        IMSIsotopeDistribution  getIsotopeDistribution() except + nogil  # wrap-doc:Gets element's isotope distribution
        void setIsotopeDistribution(IMSIsotopeDistribution & isotopes) except + nogil  # wrap-doc:Sets element's isotope distribution
        bool operator==(IMSElement & element) except + nogil 
        bool operator!=(IMSElement & element) except + nogil 

