from Types cimport *
from BaseFeature cimport *
from MSExperiment cimport *
from MSSpectrum cimport*
from SiriusMSFile cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair



cdef extern from "<OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>" namespace "OpenMS":

    cdef cppclass MetaboTargetedAssay "OpenMS::MetaboTargetedAssay":

       MetaboTargetedAssay() nogil except +
       MetaboTargetedAssay(MetaboTargetedAssay) nogil except + #wrap-ignore



       libcpp_vector[MetaboTargetedAssay] extractMetaboTargetedAssay(MSExperiment spectra,
                                                                     libcpp_map[ BaseFeature,libcpp_vector[size_t] ] feature_ms2_spectra_map,
                                                                     double precursor_rt_tol,
                                                                     double precursor_mz_distance,
                                                                     double cosine_sim_threshold,
                                                                     double transition_threshold,
                                                                     bool method_consensus_spectrum,
                                                                     bool exclude_ms2_precursor,
                                                                     unsigned int file_counter) nogil except +

       libcpp_vector[MetaboTargetedAssay] extractMetaboTargetedAssayFragmentAnnotation(libcpp_vector[ libcpp_pair[SiriusMSFile_CompoundInfo,MSSpectrum] ] v_cmp_spec,
                                                                                       double transition_threshold,
                                                                                       bool use_exact_mass,
                                                                                       bool exclude_ms2_precursor,
                                                                                       unsigned int file_counter) nogil except +