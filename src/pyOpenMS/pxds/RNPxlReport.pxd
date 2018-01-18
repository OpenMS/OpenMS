from Types cimport *

from RNPxlMarkerIonExtractor cimport *
# from ListUtilsIO cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlReport "OpenMS::RNPxlReport":
        RNPxlReport(RNPxlReport) nogil except + #wrap-ignore
        libcpp_vector[ RNPxlReportRow ] annotate(MSExperiment & spectra,
                                                 libcpp_vector[ PeptideIdentification ] & peptide_ids,
                                                 double marker_ions_tolerance) nogil except +


cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlReportRowHeader "OpenMS::RNPxlReportRowHeader":
        RNPxlReportRowHeader(RNPxlReportRowHeader) nogil except + #wrap-ignore
        String getString(String & separator) nogil except +

cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlReportRow "OpenMS::RNPxlReportRow":
        RNPxlReportRow(RNPxlReportRow) nogil except + #wrap-ignore
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
        # RNPxlMarkerIonExtractor::MarkerIonsType marker_ions
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] marker_ions
        double m_H
        double m_2H
        double m_3H
        double m_4H
        String fragment_annotation_string
        String getString(String & separator) nogil except +

