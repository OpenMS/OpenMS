from Types cimport *
from MSSpectrum cimport *
from MSExperiment cimport *
from Peak1D cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from SpectrumLookup cimport *

cdef extern from "<OpenMS/METADATA/SpectrumMetaDataLookup.h>" namespace "OpenMS":

    cdef cppclass SpectrumMetaDataLookup(SpectrumLookup):
        # wrap-inherits:
        #  SpectrumLookup

        SpectrumMetaDataLookup() nogil except +
        # SpectrumMetaDataLookup(SpectrumMetaDataLookup) nogil except + # private

        void readSpectra(MSExperiment spectra, String scan_regexp, bool get_precursor_rt) nogil except +

        void getSpectrumMetaData(Size index, SpectrumMetaData& meta) nogil except +

        void getSpectrumMetaData(String spectrum_ref, SpectrumMetaData& meta) nogil except +

        void getSpectrumMetaData(String spectrum_ref, SpectrumMetaData& meta, unsigned char flags) nogil except +

        void setSpectraDataRef(const String & spectra_data) nogil except +

#
## wrap static methods
#
cdef extern from "<OpenMS/METADATA/SpectrumMetaDataLookup.h>" namespace "OpenMS::SpectrumMetaDataLookup":

    ###   TODO: Missing optional parameters:
    ###   boost::regex& scan_regexp 
    ###   std::map<Size, double>& precursor_rts 
    void getSpectrumMetaData(MSSpectrum spectrum,
                             SpectrumMetaData& meta) nogil except + # wrap-attach:SpectrumMetaDataLookup

    bool addMissingRTsToPeptideIDs(libcpp_vector[PeptideIdentification], 
                                   String filename, bool stop_on_error) nogil except + # wrap-attach:SpectrumMetaDataLookup

    bool addMissingSpectrumReferences(libcpp_vector[PeptideIdentification], 
                                   String filename, bool stop_on_error, 
                                   bool override_spectra_data, 
                                   bool override_spectra_references,
                                   libcpp_vector[ProteinIdentification] proteins) nogil except + # wrap-attach:SpectrumMetaDataLookup

cdef extern from "<OpenMS/METADATA/SpectrumMetaDataLookup.h>" namespace "OpenMS::SpectrumMetaDataLookup":

    cdef cppclass SpectrumMetaData "OpenMS::SpectrumMetaDataLookup::SpectrumMetaData":

      SpectrumMetaData() nogil except +
      SpectrumMetaData(SpectrumMetaData) nogil except + # wrap-ignore

      double rt
      double precursor_rt
      double precursor_mz
      Int precursor_charge
      Size ms_level
      Int scan_number
      String native_id

