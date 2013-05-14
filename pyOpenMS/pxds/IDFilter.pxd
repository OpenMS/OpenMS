from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool

from String cimport *
from DefaultParamHandler cimport *

from PeptideIdentification cimport *
from ProteinIdentification cimport *
from FASTAFile cimport *

from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/FILTERING/ID/IDFilter.h>" namespace "OpenMS":

    cdef cppclass IDFilter:

        IDFilter()           nogil except +
        IDFilter(IDFilter)   nogil except + # wrap-ignore

        void filterIdentificationsByThreshold(PeptideIdentification& identification, DoubleReal threshold_fraction, PeptideIdentification& filtered_identification)
        void filterIdentificationsByScore(PeptideIdentification& identification, DoubleReal threshold_score, PeptideIdentification& filtered_identification)
        void filterIdentificationsByBestNHits(PeptideIdentification& identification, Size n, PeptideIdentification& filtered_identification)
        void filterIdentificationsByBestNToMHits(PeptideIdentification& identification, Size n, Size m, PeptideIdentification& filtered_identification)

        void filterIdentificationsByThreshold(ProteinIdentification& identification, DoubleReal threshold_fraction, ProteinIdentification& filtered_identification)
        void filterIdentificationsByScore(ProteinIdentification& identification, DoubleReal threshold_score, ProteinIdentification& filtered_identification)
        void filterIdentificationsByBestNHits(ProteinIdentification& identification, Size n, ProteinIdentification& filtered_identification)
        void filterIdentificationsByBestNToMHits(ProteinIdentification& identification, Size n, Size m, ProteinIdentification& filtered_identification)

        void filterIdentificationsByBestHits(PeptideIdentification& identification, PeptideIdentification& filtered_identification, bool strict)

        void filterIdentificationsByProteins(PeptideIdentification& identification, libcpp_vector[FASTAEntry]& proteins, PeptideIdentification& filtered_identification, bool no_protein_identifiers)
        void filterIdentificationsByProteins(ProteinIdentification& identification, libcpp_vector[FASTAEntry]& proteins, ProteinIdentification& filtered_identification)
        void filterIdentificationsByExclusionPeptides(PeptideIdentification& identification, libcpp_set[String]& peptides, PeptideIdentification& filtered_identification)
        void filterIdentificationsByLength(PeptideIdentification& identification, PeptideIdentification& filtered_identification, Size min_length, Size max_length)
        void filterIdentificationsByCharge(PeptideIdentification& identification, Int charge, PeptideIdentification& filtered_identification)

        void filterIdentificationsByVariableModifications(PeptideIdentification& identification, libcpp_vector[String]& fixed_modifications, PeptideIdentification& filtered_identification)
        void removeUnreferencedProteinHits(ProteinIdentification& identification, libcpp_vector[PeptideIdentification] peptide_identifications, ProteinIdentification& filtered_identification)
        void filterIdentificationsUnique(PeptideIdentification& identification, PeptideIdentification& filtered_identification)
        void filterIdentificationsByMzError(PeptideIdentification& identification, DoubleReal mass_error, bool unit_ppm, PeptideIdentification& filtered_identification)
        void filterIdentificationsByRTPValues(PeptideIdentification& identification, PeptideIdentification& filtered_identification, DoubleReal p_value)
        void filterIdentificationsByRTFirstDimPValues(PeptideIdentification& identification, PeptideIdentification& filtered_identification, DoubleReal p_value) 
        void filterIdentificationsByThresholds(MSExperiment[Peak1D,ChromatogramPeak]& experiment, DoubleReal peptide_threshold_fraction, DoubleReal protein_threshold_fraction)
        void filterIdentificationsByScores(MSExperiment[Peak1D,ChromatogramPeak]& experiment, DoubleReal peptide_threshold_score, DoubleReal protein_threshold_score)
        void filterIdentificationsByBestNHits(MSExperiment[Peak1D,ChromatogramPeak]& experiment, Size n)
        void filterIdentificationsByProteins(MSExperiment[Peak1D,ChromatogramPeak]& experiment, libcpp_vector[FASTAEntry]& proteins)

