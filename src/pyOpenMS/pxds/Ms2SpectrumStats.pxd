from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from String cimport *
from QCBase cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/QC/Ms2SpectrumStats.h>" namespace "OpenMS":
    
    cdef cppclass Ms2SpectrumStats(QCBase):
        # wrap-inherits:
        #    QCBase
        #
        # wrap-doc:
        #  QC metric to determine the number of MS2 scans per MS1 scan over RT
        #  
        #  Ms2SpectrumStats collects data from MS2 scans and stores the result into 
        #  PeptideIdentifications, which already exist in the FeatureMap, or are newly 
        #  created as empty PeptideIdentifications (with no sequence).

        Ms2SpectrumStats() except + nogil 
        Ms2SpectrumStats(Ms2SpectrumStats &) except + nogil 

        PeptideIdentificationList compute(MSExperiment & exp, FeatureMap & features, SpectraMap & map_to_spectrum) except + nogil 
            # wrap-doc:
            #  Calculate the ScanEventNumber, find all unidentified MS2-Spectra and add 
            #  them to unassigned PeptideIdentifications, write meta values "ScanEventNumber" 
            #  and "identified" in PeptideIdentification.
            #  
            #  :param exp: Imported calibrated MzML file as MSExperiment
            #  :param features: Imported featureXML file after FDR as FeatureMap
            #  :param map_to_spectrum: Map to find index of spectrum given by meta value at PepID
            #  :return: Unassigned peptide identifications newly generated from unidentified MS2-Spectra

        const String & getName() except + nogil  # wrap-doc:Returns the name of the metric

cdef extern from "<OpenMS/QC/Ms2SpectrumStats.h>" namespace "OpenMS::Ms2SpectrumStats":
    
    cdef cppclass ScanEvent "OpenMS::Ms2SpectrumStats::ScanEvent":
        # wrap-doc:
        #  Structure for storing scan event information
        
        ScanEvent(UInt32 sem, bool ms2) except + nogil 
        ScanEvent(ScanEvent &) except + nogil 
        
        UInt32 scan_event_number  # wrap-doc:Scan event number
        bool ms2_presence  # wrap-doc:MS2 presence
