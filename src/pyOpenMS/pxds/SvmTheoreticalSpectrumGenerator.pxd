from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from SimTypes cimport *
from SVMWrapper cimport *
from AASequence cimport *

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>" namespace "OpenMS":
    
    cdef cppclass SvmTheoreticalSpectrumGenerator "OpenMS::SvmTheoreticalSpectrumGenerator":
        SvmTheoreticalSpectrumGenerator() except + nogil 
        SvmTheoreticalSpectrumGenerator(SvmTheoreticalSpectrumGenerator &) except + nogil 
        # void simulate(MSSpectrum &spectrum, AASequence &peptide, boost::random::mt19937_64& rng, Size precursor_charge) except + nogil 
        void load() except + nogil 
        libcpp_vector[ IonType ]  getIonTypes() except + nogil 

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>" namespace "OpenMS::SvmTheoreticalSpectrumGenerator":
    
    cdef cppclass SvmModelParameterSet "OpenMS::SvmTheoreticalSpectrumGenerator::SvmModelParameterSet":

        SvmModelParameterSet(SvmModelParameterSet) except + nogil  #wrap-ignore

        # libcpp_vector[ shared_ptr[ SVMWrapper ] ] class_models
        # libcpp_vector[ shared_ptr[ SVMWrapper ] ] reg_models
        # TODO STL map with wrapped key
        # libcpp_map[ ResidueType, double ] _intensities
        # libcpp_vector[ IonType ] ion_types
        # TODO STL map with wrapped key
        # libcpp_map[ IonType, libcpp_vector[ IonType ] ] secondary_types
        Size number_intensity_levels
        Size number_regions
        libcpp_vector[ double ] feature_max
        libcpp_vector[ double ] feature_min
        double scaling_lower
        double scaling_upper
        libcpp_vector[ double ] intensity_bin_boarders
        libcpp_vector[ double ] intensity_bin_values
        # TODO nested STL
        # libcpp_map[ libcpp_pair[ IonType, Size ], libcpp_vector[ libcpp_vector[ double ] ] ] conditional_prob

cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>" namespace "OpenMS::SvmTheoreticalSpectrumGenerator":
    
    cdef cppclass IonType "OpenMS::SvmTheoreticalSpectrumGenerator::IonType":
        IonType() except + nogil  # compiler
        IonType(IonType &) except + nogil  # compiler
        IonType(ResidueType residue, EmpiricalFormula l, Int charge) except + nogil 
        bool operator<(IonType &rhs) except + nogil 
        ResidueType residue
        EmpiricalFormula loss
        Int charge

# cdef extern from "<OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGenerator.h>" namespace "OpenMS::SvmTheoreticalSpectrumGenerator":
#    
#    cdef cppclass DescriptorSet "OpenMS::SvmTheoreticalSpectrumGenerator::DescriptorSet":
#        DescriptorSet(DescriptorSet) except + nogil  #wrap-ignore
#        libcpp_vector[svm_node] descriptors # we would have to wrap libsvm for this
