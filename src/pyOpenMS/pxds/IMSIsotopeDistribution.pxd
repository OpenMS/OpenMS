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
        # wrap-doc:
                #  Represents a distribution of isotopes restricted to the first K elements
                #  
                #  Represents a distribution of isotopes of chemical elements as a list
                #  of peaks each as a pair of mass and abundance. 'IsotopeDistribution'
                #  unlike 'IsotopeSpecies' has one abundance per a nominal mass.
                #  Here is an example in the format (mass; abundance %)
                #  for molecule H2O (values are taken randomly):
                #  
                #  - IsotopeDistribution
                #      (18.00221; 99.03 %)
                #      (19.00334; 0.8 %)
                #      (20.00476; 0.17 %)
                #  
                #  - IsotopeSpecies
                #      (18.00197; 98.012 %)
                #      (18.00989; 1.018 %)
                #      (19.00312; 0.683 %)
                #      (19.00531; 0.117 %)
                #      (20.00413; 0.134 %)
                #      (20.00831; 0.036 %)
                #  
                #  To the sake of faster computations distribution is restricted
                #  to the first K elements, where K can be set by adjusting size
                #  'SIZE' of distribution. @note For the elements most abundant in
                #  living beings (CHNOPS) this restriction is negligible, since abundances
                #  decrease dramatically in isotopes order and are usually of no interest
                #  starting from +10 isotope.
                #  
                #  'IsotopeDistribution' implements folding with other distribution using an
                #  algorithm described in details in paper:
                #  Boecker et al. "Decomposing metabolic isotope patterns" WABI 2006. doi: 10.1007/11851561_2
                #  
                #  Folding with itself is done using Russian Multiplication Scheme
                
        IMSIsotopeDistribution() except + nogil 
        IMSIsotopeDistribution(IMSIsotopeDistribution &) except + nogil 

        abundance_type ABUNDANCES_SUM_ERROR
        size_type SIZE
        IMSIsotopeDistribution(nominal_mass_type nominalMass) except + nogil 
        IMSIsotopeDistribution(mass_type mass) except + nogil 
        IMSIsotopeDistribution(libcpp_vector[IMSIsotopeDistribution_Peak] & peaks, nominal_mass_type nominalMass) except + nogil 
        size_type size() except + nogil 
        bool operator==(IMSIsotopeDistribution & distribution) except + nogil 
        bool operator!=(IMSIsotopeDistribution & distribution) except + nogil 
        # IMSIsotopeDistribution operator*=(IMSIsotopeDistribution & distribution) except + nogil 
        # IMSIsotopeDistribution operator*=(unsigned int pow_) except + nogil 
        mass_type getMass(size_type i) except + nogil 
        abundance_type getAbundance(size_type i) except + nogil 
        mass_type getAverageMass() except + nogil 
        nominal_mass_type getNominalMass() except + nogil 
        void setNominalMass(nominal_mass_type nominalMass) except + nogil 
        masses_container getMasses() except + nogil  # wrap-doc:Gets a mass of isotope 'i'
        libcpp_vector[abundance_type] getAbundances() except + nogil  # wrap-doc:Gets an abundance of isotope 'i'
        void normalize() except + nogil  # wrap-doc:Normalizes distribution, i.e. scaling abundances to be summed up to 1 with an error
        bool empty() except + nogil  # wrap-doc:Returns true if the distribution has no peaks, false - otherwise

cdef extern from "<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>" namespace "OpenMS::ims::IMSIsotopeDistribution":
    
    cdef cppclass IMSIsotopeDistribution_Peak "OpenMS::ims::IMSIsotopeDistribution::Peak":
        IMSIsotopeDistribution_Peak(IMSIsotopeDistribution_Peak) except + nogil  #wrap-ignore
        mass_type mass
        abundance_type abundance
        IMSIsotopeDistribution_Peak(mass_type mass, abundance_type abundance) except + nogil 
        bool operator==(IMSIsotopeDistribution_Peak & peak) except + nogil 

