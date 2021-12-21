from Types cimport *
from SiriusMSFile cimport *
from SiriusFragmentAnnotation cimport *
from FeatureMapping cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.unordered_map cimport unordered_map as libcpp_unordered_map

cdef extern from "<OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>" namespace "OpenMS":

    cdef cppclass MetaboTargetedAssay "OpenMS::MetaboTargetedAssay":

       MetaboTargetedAssay() nogil except + # wrap-doc:This class provides methods for the extraction of targeted assays for metabolomics
       MetaboTargetedAssay(MetaboTargetedAssay &) nogil except + # compiler

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
            # wrap-doc:
                #   Extract a vector of MetaboTargetedAssays without using fragment annotation
                #   -----
                #   :param spectra: Input of MSExperiment with spectra information
                #   :param feature_ms2_spectra_map: FeatureMapping class with associated MS2 spectra
                #   :param precursor_rt_tol: Retention time tolerance of the precursor
                #   :param precursor_mz_distance: Max m/z distance of the precursor entries of two spectra to be merged
                #   :param cosine_sim_threshold: Cosine similarty threshold for the usage of SpectraMerger
                #   :param transition_threshold: Intensity threshold for MS2 peak used in MetaboTargetedAssay
                #   :param min_fragment_mz: Minimum m/z a fragment ion has to have to be considered as a transition
                #   :param max_fragment_mz: Maximum m/z a fragment ion has to have to be considered as a transition
                #   :param method_consensus_spectrum: Boolean to use consensus spectrum method
                #   :param exclude_ms2_precursor: Boolean to exclude MS2 precursor from MetaboTargetedAssay
                #   :param file_counter: Count if multiple files are used
                #   :returns: Vector of MetaboTargetedAssay

       libcpp_vector[ MetaboTargetedAssay ] extractMetaboTargetedAssayFragmentAnnotation(libcpp_vector[ MetaboTargetedAssay_CompoundTargetDecoyPair ]& v_cmp_spec,
                                                                                         double& transition_threshold,
                                                                                         double& min_fragment_mz,
                                                                                         double& max_fragment_mz,
                                                                                         bool& use_exact_mass,
                                                                                         bool& exclude_ms2_precursor,
                                                                                         unsigned int& file_counter) nogil except +
            # wrap-doc:
                #   Extract a vector of MetaboTargetedAssays using fragment annotation
                #   -----
                #   :param v_cmp_spec: Vector of CompoundInfo with associated fragment annotated MSspectrum
                #   :param transition_threshold: Intensity threshold for MS2 peak used in MetaboTargetedAssay
                #   :param min_fragment_mz: Minimum m/z a fragment ion has to have to be considered as a transition
                #   :param max_fragment_mz: Maximum m/z a fragment ion has to have to be considered as a transition
                #   :param use_exact_mass: Boolean if exact mass should be used as peak mass for annotated fragments
                #   :param exclude_ms2_precursor: Boolean to exclude MS2 precursor from MetaboTargetedAssay
                #   :param file_counter: Count if multiple files are used.
                #   :returns: Vector of MetaboTargetedAssay

       libcpp_vector[ MetaboTargetedAssay_CompoundTargetDecoyPair ] pairCompoundWithAnnotatedSpectra(libcpp_vector[ SiriusMSFile_CompoundInfo ]& v_cmpinfo,
                                                                                                     libcpp_vector[ SiriusFragmentAnnotation_SiriusTargetDecoySpectra]& annotated_spectra) nogil except +
            # wrap-doc:
                #   Pair compound information (SiriusMSFile) with the annotated target and decoy spectrum from SIRIUS/Passatutto based on the m_id (unique identifier composed of description_filepath_native_id_k introduced in the SiriusMSConverter)
                #   -----
                #   :param v_cmpinfo: Vector of SiriusMSFile::CompoundInfo
                #   :param annotated_spectra: Vector of SiriusTargetDecoySpectra
                #   :return: Vector of MetaboTargetedAssay::CompoundTargetDecoyPair

#       libcpp_unordered_map[ UInt64, libcpp_vector[ MetaboTargetedAssay ] ] buildAmbiguityGroup(libcpp_vector[ MetaboTargetedAssay]& v_mta,
#                                                                                              double& ar_mz_tol,
#                                                                                              double& ar_rt_tol,
#                                                                                              String& ar_mz_tol_unit_res,
#                                                                                              size_t in_files_size) nogil except +

#        void resolveAmbiguityGroup(libcpp_unordered_map[ UInt64, libcpp_vector[ MetaboTargetedAssay ] ]& map_mta_filter,
#                                   double& total_occurrence_filter,
#                                   size_t in_files_size) nogil except +

    cdef cppclass MetaboTargetedAssay_CompoundTargetDecoyPair "OpenMS::MetaboTargetedAssay::CompoundTargetDecoyPair":

        MetaboTargetedAssay_CompoundTargetDecoyPair() nogil except +
        MetaboTargetedAssay_CompoundTargetDecoyPair(MetaboTargetedAssay_CompoundTargetDecoyPair &) nogil except + # compiler
