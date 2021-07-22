from Types cimport *

from NuXLMarkerIonExtractor cimport *
# from ListUtilsIO cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/NUXL/NuXLReport.h>" namespace "OpenMS":
    
    cdef cppclass NuXLReport "OpenMS::NuXLReport":
        NuXLReport(NuXLReport) nogil except + #wrap-ignore
        libcpp_vector[ NuXLReportRow ] annotate(MSExperiment & spectra,
                                                 libcpp_vector[ PeptideIdentification ] & peptide_ids,
                                                 StringList meta_values_to_export,
                                                 double marker_ions_tolerance) nogil except +


cdef extern from "<OpenMS/ANALYSIS/NUXL/NuXLReport.h>" namespace "OpenMS":
    
    cdef cppclass NuXLReportRowHeader "OpenMS::NuXLReportRowHeader":
        NuXLReportRowHeader(NuXLReportRowHeader) nogil except + #wrap-ignore
        String getString(const String & separator, StringList meta_values_to_export) nogil except +

cdef extern from "<OpenMS/ANALYSIS/NUXL/NuXLReport.h>" namespace "OpenMS":
    
    cdef cppclass NuXLReportRow "OpenMS::NuXLReportRow":
        NuXLReportRow(NuXLReportRow) nogil except + #wrap-ignore

        bool no_id
        double rt
        double original_mz
        String accessions
        String peptide
        String NA
        Int charge
        double score
        int rank
        double best_localization_score
        String localization_scores
        String best_localization
        double peptide_weight
        double NA_weight
        double xl_weight
        StringList meta_values
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] marker_ions # NuXLMarkerIonExtractor::MarkerIonsType marker_ions
        double abs_prec_error
        double rel_prec_error
        double m_H
        double m_2H
        double m_3H
        double m_4H
        String fragment_annotation
        String getString(const String & separator) nogil except +

