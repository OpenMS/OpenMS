from Types cimport *
from libcpp cimport bool
from IMSIsotopeDistribution cimport *
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>" namespace "OpenMS::ims::IMSElement":

    ctypedef libcpp_string name_type
    ctypedef IMSIsotopeDistribution isotopes_type

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>" namespace "OpenMS::ims":
    
    cdef cppclass IMSElement "OpenMS::ims::IMSElement":
        IMSElement() nogil except +
        IMSElement(IMSElement) nogil except +
        # mass_type ELECTRON_MASS_IN_U read-only
        IMSElement(name_type & name, isotopes_type & isotopes) nogil except +
        IMSElement(name_type & name, mass_type mass) nogil except +
        IMSElement(name_type & name, nominal_mass_type nominal_mass) nogil except +
        name_type  getName() nogil except +
        void setName(name_type & name) nogil except +
        name_type  getSequence() nogil except +
        void setSequence(name_type & sequence) nogil except +
        nominal_mass_type getNominalMass() nogil except +
        mass_type getMass(size_type index) nogil except +
        mass_type getAverageMass() nogil except +
        mass_type getIonMass(int electrons_number) nogil except +
        IMSIsotopeDistribution  getIsotopeDistribution() nogil except +
        void setIsotopeDistribution(IMSIsotopeDistribution & isotopes) nogil except +
        bool operator==(IMSElement & element) nogil except +
        bool operator!=(IMSElement & element) nogil except +

