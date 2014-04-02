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

        void filterIdentificationsByThreshold(PeptideIdentification& identification, double threshold_fraction, PeptideIdentification& filtered_identification) nogil except +
        void filterIdentificationsByScore(PeptideIdentification& identification, double threshold_score, PeptideIdentification& filtered_identification) nogil except +
        void filterIdentificationsByBestNHits(PeptideIdentification& identification, Size n, PeptideIdentification& filtered_identification) nogil except +
        void filterIdentificationsByBestNToMHits(PeptideIdentification& identification, Size n, Size m, PeptideIdentification& filtered_identification) nogil except +

        void filterIdentificationsByThreshold(ProteinIdentification& identification, double threshold_fraction, ProteinIdentification& filtered_identification) nogil except +
        void filterIdentificationsByScore(ProteinIdentification& identification, double threshold_score, ProteinIdentification& filtered_identification) nogil except +
        void filterIdentificationsByBestNHits(ProteinIdentification& identification, Size n, ProteinIdentification& filtered_identification) nogil except +
        void filterIdentificationsByBestNToMHits(ProteinIdentification& identification, Size n, Size m, ProteinIdentification& filtered_identification) nogil except +

        void filterIdentificationsByBestHits(PeptideIdentification& identification, PeptideIdentification& filtered_identification, bool strict) nogil except +

        void filterIdentificationsByProteins(PeptideIdentification& identification, libcpp_vector[FASTAEntry]& proteins, PeptideIdentification& filtered_identification, bool no_protein_identifiers) nogil except +
        void filterIdentificationsByProteins(ProteinIdentification& identification, libcpp_vector[FASTAEntry]& proteins, ProteinIdentification& filtered_identification) nogil except +
        void filterIdentificationsByExclusionPeptides(PeptideIdentification& identification, libcpp_set[String]& peptides, PeptideIdentification& filtered_identification) nogil except +
        void filterIdentificationsByLength(PeptideIdentification& identification, PeptideIdentification& filtered_identification, Size min_length, Size max_length) nogil except +
        void filterIdentificationsByCharge(PeptideIdentification& identification, Int charge, PeptideIdentification& filtered_identification) nogil except +

        void filterIdentificationsByVariableModifications(PeptideIdentification& identification, libcpp_vector[String]& fixed_modifications, PeptideIdentification& filtered_identification) nogil except +
        void removeUnreferencedProteinHits(ProteinIdentification& identification, libcpp_vector[PeptideIdentification] peptide_identifications, ProteinIdentification& filtered_identification) nogil except +
        void filterIdentificationsUnique(PeptideIdentification& identification, PeptideIdentification& filtered_identification) nogil except +
        void filterIdentificationsByMzError(PeptideIdentification& identification, double mass_error, bool unit_ppm, PeptideIdentification& filtered_identification) nogil except +
        void filterIdentificationsByRTPValues(PeptideIdentification& identification, PeptideIdentification& filtered_identification, double p_value) nogil except +
        void filterIdentificationsByRTFirstDimPValues(PeptideIdentification& identification, PeptideIdentification& filtered_identification, double p_value)  nogil except +
        void filterIdentificationsByThresholds(MSExperiment[Peak1D,ChromatogramPeak]& experiment, double peptide_threshold_fraction, double protein_threshold_fraction) nogil except +
        void filterIdentificationsByScores(MSExperiment[Peak1D,ChromatogramPeak]& experiment, double peptide_threshold_score, double protein_threshold_score) nogil except +
        void filterIdentificationsByBestNHits(MSExperiment[Peak1D,ChromatogramPeak]& experiment, Size n) nogil except +
        void filterIdentificationsByProteins(MSExperiment[Peak1D,ChromatogramPeak]& experiment, libcpp_vector[FASTAEntry]& proteins) nogil except +

