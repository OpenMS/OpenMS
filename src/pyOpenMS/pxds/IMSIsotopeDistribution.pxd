from Types cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>" namespace "OpenMS::ims::IMSIsotopeDistribution":

    ctypedef signed int size_type
    ctypedef double mass_type
    ctypedef double abundance_type
    ctypedef unsigned int nominal_mass_type
    ctypedef libcpp_vector[mass_type] masses_container
    ctypedef IMSIsotopeDistribution_Peak peak_type
    ctypedef libcpp_vector[peak_type] peaks_container
    ctypedef libcpp_vector[abundance_type] abundances_container
    
cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>" namespace "OpenMS::ims":

    cdef cppclass IMSIsotopeDistribution "OpenMS::ims::IMSIsotopeDistribution":
        IMSIsotopeDistribution() nogil except +
        IMSIsotopeDistribution(IMSIsotopeDistribution) nogil except +
        abundance_type ABUNDANCES_SUM_ERROR
        size_type SIZE
        IMSIsotopeDistribution(nominal_mass_type nominalMass) nogil except +
        IMSIsotopeDistribution(mass_type mass) nogil except +
        IMSIsotopeDistribution(libcpp_vector[IMSIsotopeDistribution_Peak] & peaks, nominal_mass_type nominalMass) nogil except +
        size_type size() nogil except +
        bool operator==(IMSIsotopeDistribution & distribution) nogil except +
        bool operator!=(IMSIsotopeDistribution & distribution) nogil except +
        # IMSIsotopeDistribution operator*=(IMSIsotopeDistribution & distribution) nogil except +
        # IMSIsotopeDistribution operator*=(unsigned int pow_) nogil except +
        mass_type getMass(size_type i) nogil except +
        abundance_type getAbundance(size_type i) nogil except +
        mass_type getAverageMass() nogil except +
        nominal_mass_type getNominalMass() nogil except +
        void setNominalMass(nominal_mass_type nominalMass) nogil except +
        masses_container getMasses() nogil except + 
        libcpp_vector[abundance_type] getAbundances() nogil except +
        void normalize() nogil except +
        bool empty() nogil except +

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>" namespace "OpenMS::ims::IMSIsotopeDistribution":
    
    cdef cppclass IMSIsotopeDistribution_Peak "OpenMS::ims::IMSIsotopeDistribution::Peak":
        IMSIsotopeDistribution_Peak(IMSIsotopeDistribution_Peak) nogil except + #wrap-ignore
        mass_type mass
        abundance_type abundance
        IMSIsotopeDistribution_Peak(mass_type mass, abundance_type abundance) nogil except +
        bool operator==(IMSIsotopeDistribution_Peak & peak) nogil except +

