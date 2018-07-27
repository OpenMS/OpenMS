from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool
from XLPrecursor cimport *
from AASeqWithMass cimport *
from DoubleList cimport *
from StringList cimport *
from ResidueModification cimport *
from FASTAFile cimport *
from EnzymaticDigestion cimport *
from ProteinProteinCrossLink cimport *
from CrossLinkSpectrumMatch cimport *
from PeptideHit cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from MSExperiment cimport *
from EnzymaticDigestion cimport *


cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLHelper.h>" namespace "OpenMS":

    cdef cppclass OPXLHelper:

        OPXLHelper(OPXLHelper) nogil except +
        OPXLHelper() nogil except +

        libcpp_vector[ XLPrecursor ] enumerateCrossLinksAndMasses(libcpp_vector[ AASeqWithMass ]  peptides,
                                                                  double cross_link_mass_light,
                                                                  DoubleList cross_link_mass_mono_link,
                                                                  StringList cross_link_residue1,
                                                                  StringList cross_link_residue2,
                                                                  libcpp_vector[ double ]& spectrum_precursors,
                                                                  libcpp_vector[ int ]& precursor_correction_positions,
                                                                  double precursor_mass_tolerance,
                                                                  bool precursor_mass_tolerance_unit_ppm)  nogil except +

        libcpp_vector[ ResidueModification ] getModificationsFromStringList(StringList modNames)  nogil except +

        libcpp_vector[ AASeqWithMass ] digestDatabase(libcpp_vector[ FASTAEntry ] fasta_db,
                                                      EnzymaticDigestion digestor,
                                                      Size min_peptide_length,
                                                      StringList cross_link_residue1,
                                                      StringList cross_link_residue2,
                                                      libcpp_vector[ ResidueModification ] fixed_modifications,
                                                      libcpp_vector[ ResidueModification ] variable_modifications,
                                                      Size max_variable_mods_per_peptide) nogil except +

        libcpp_vector[ ProteinProteinCrossLink ] buildCandidates(libcpp_vector[ XLPrecursor ]& candidates,
                                                                 libcpp_vector[ int ]& precursor_corrections,
                                                                 libcpp_vector[ int ]& precursor_correction_positions,
                                                                 libcpp_vector[ AASeqWithMass ]& peptide_masses,
                                                                 const StringList& cross_link_residue1,
                                                                 const StringList& cross_link_residue2,
                                                                 double cross_link_mass,
                                                                 DoubleList cross_link_mass_mono_link,
                                                                 libcpp_vector[ double ]& spectrum_precursor_vector,
                                                                 libcpp_vector[ double ]& allowed_error_vector,
                                                                 String cross_link_name) nogil except +


        void buildFragmentAnnotations(libcpp_vector[ PeptideHit_PeakAnnotation ]& frag_annotations,
                                      libcpp_vector[ libcpp_pair[ size_t, size_t ] ] matching,
                                      MSSpectrum theoretical_spectrum,
                                      MSSpectrum experiment_spectrum) nogil except +

        void buildPeptideIDs(libcpp_vector[ PeptideIdentification ]& peptide_ids,
                             libcpp_vector[ CrossLinkSpectrumMatch ] top_csms_spectrum,
                             libcpp_vector[ libcpp_vector[ CrossLinkSpectrumMatch ] ]& all_top_csms,
                             Size all_top_csms_current_index,
                             MSExperiment spectra,
                             Size scan_index,
                             Size scan_index_heavy) nogil except +


        void addProteinPositionMetaValues(libcpp_vector[ PeptideIdentification ]& peptide_ids) nogil except +
