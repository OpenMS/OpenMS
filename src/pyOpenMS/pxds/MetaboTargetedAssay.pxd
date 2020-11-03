from Types cimport *
from BaseFeature cimport *
from MSExperiment cimport *
from MSSpectrum cimport *
from SiriusMSFile cimport *
from SiriusFragmentAnnotation cimport *
from FeatureMapping cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair

cdef extern from "<OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>" namespace "OpenMS":

    cdef cppclass MetaboTargetedAssay "OpenMS::MetaboTargetedAssay":

       MetaboTargetedAssay() nogil except +
       MetaboTargetedAssay(MetaboTargetedAssay) nogil except +

       libcpp_vector[ MetaboTargetedAssay ] extractMetaboTargetedAssay(MSExperiment& spectra,
                                                                       FeatureMapping_FeatureToMs2Indices& feature_ms2_index,
                                                                       double& precursor_rt_tol,
                                                                       double& precursor_mz_distance,
                                                                       double& cosine_sim_threshold,
                                                                       double& transition_threshold,
                                                                       double& min_fragment_mz,
                                                                       double& max_fragment_mz,
                                                                       bool& method_consensus_spectrum,
                                                                       bool& exclude_ms2_precursor,
                                                                       unsigned int& file_counter) nogil except +

       libcpp_vector[ MetaboTargetedAssay ] extractMetaboTargetedAssayFragmentAnnotation(libcpp_vector[ MetaboTargetedAssay_CompoundTargetDecoyPair ]& v_cmp_spec,
                                                                                         double& transition_threshold,
                                                                                         double& min_fragment_mz,
                                                                                         double& max_fragment_mz,
                                                                                         bool& use_exact_mass,
                                                                                         bool& exclude_ms2_precursor,
                                                                                         unsigned int& file_counter) nogil except +
       libcpp_vector[ MetaboTargetedAssay_CompoundTargetDecoyPair ] pairCompoundWithAnnotatedSpectra(libcpp_vector[ SiriusMSFile_CompoundInfo ]& v_cmpinfo, libcpp_vector[ SiriusFragmentAnnotation_SiriusTargetDecoySpectra]& annotated_spectra) nogil except +

    cdef cppclass MetaboTargetedAssay_CompoundTargetDecoyPair "OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair":

        MetaboTargetedAssay_CompoundTargetDecoyPair() nogil except +
        MetaboTargetedAssay_CompoundTargetDecoyPair(MetaboTargetedAssay_CompoundTargetDecoyPair) nogil except +
