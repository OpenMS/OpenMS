from Types cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>" namespace "OpenMS":

    cdef cppclass MRMFQC_ComponentGroupQCs "OpenMS::MRMFeatureQC::ComponentGroupQCs":
        MRMFQC_ComponentGroupQCs()
        MRMFQC_ComponentGroupQCs(MRMFQC_ComponentGroupQCs &) # no-wrap

        String component_group_name
        double retention_time_l
        double retention_time_u
        double intensity_l
        double intensity_u
        double overall_quality_l
        double overall_quality_u
        Int n_heavy_l
        Int n_heavy_u
        Int n_light_l
        Int n_light_u
        Int n_detecting_l
        Int n_detecting_u
        Int n_quantifying_l
        Int n_quantifying_u
        Int n_identifying_l
        Int n_identifying_u
        Int n_transitions_l
        Int n_transitions_u
        String ion_ratio_pair_name_1
        String ion_ratio_pair_name_2
        double ion_ratio_l
        double ion_ratio_u
        String ion_ratio_feature_name
        # libcpp_map[String,libcpp_pair[double,double]] meta_value_qc
        ## currently not supported

    cdef cppclass MRMFQC_ComponentQCs "OpenMS::MRMFeatureQC::ComponentQCs":
        MRMFQC_ComponentQCs()
        MRMFQC_ComponentQCs(MRMFQC_ComponentQCs &) # no-wrap

        String component_name 
        double retention_time_l
        double retention_time_u
        double intensity_l
        double intensity_u
        double overall_quality_l
        double overall_quality_u
        # libcpp_map[String,libcpp_pair[double,double]] meta_value_qc
        ## currently not supported

    cdef cppclass MRMFQC_ComponentGroupPairQCs "OpenMS::MRMFeatureQC::ComponentGroupPairQCs":
        MRMFQC_ComponentGroupPairQCs()
        MRMFQC_ComponentGroupPairQCs(MRMFQC_ComponentGroupPairQCs &) # no-wrap

        String component_group_name 
        String resolution_pair_name 
        double resolution_l
        double resolution_u
        double rt_diff_l
        double rt_diff_u

    cdef cppclass MRMFeatureQC:

        MRMFeatureQC() nogil except +
        MRMFeatureQC(MRMFeatureQC &) nogil except +
        
        libcpp_vector[MRMFQC_ComponentQCs] component_qcs
        libcpp_vector[MRMFQC_ComponentGroupQCs] component_group_qcs
        libcpp_vector[MRMFQC_ComponentGroupPairQCs] component_group_pair_qcs

