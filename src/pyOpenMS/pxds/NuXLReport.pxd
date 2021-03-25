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
                                                 double marker_ions_tolerance) nogil except +


cdef extern from "<OpenMS/ANALYSIS/NUXL/NuXLReport.h>" namespace "OpenMS":
    
    cdef cppclass NuXLReportRowHeader "OpenMS::NuXLReportRowHeader":
        NuXLReportRowHeader(NuXLReportRowHeader) nogil except + #wrap-ignore
        String getString(const String & separator) nogil except +

cdef extern from "<OpenMS/ANALYSIS/NUXL/NuXLReport.h>" namespace "OpenMS":
    
    cdef cppclass NuXLReportRow "OpenMS::NuXLReportRow":
        NuXLReportRow(NuXLReportRow) nogil except + #wrap-ignore

        bool no_id
        double rt
        double original_mz
        String accessions
        String RNA
        String peptide
        double best_localization_score
        String localization_scores
        String best_localization
        Int charge
        double score
        double peptide_weight
        double RNA_weight
        double xl_weight
        double abs_prec_error
        double rel_prec_error
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] marker_ions # NuXLMarkerIonExtractor::MarkerIonsType marker_ions
        double m_H
        double m_2H
        double m_3H
        double m_4H
        int rank

        String getString(const String & separator) nogil except +

