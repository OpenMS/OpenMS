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

        SpectrumMetaDataLookup() except + nogil 
        # private
        SpectrumMetaDataLookup(SpectrumMetaDataLookup) except + nogil  # wrap-ignore

        void readSpectra(MSExperiment spectra, String scan_regexp, bool get_precursor_rt) except + nogil 
        # wrap-doc:
                #  Read spectra and store their meta data
                #  
                #  :param SpectrumContainer: Spectrum container class, must support `size` and `operator[]`
                #  :param spectra: Container of spectra
                #  :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>")
                #  :param get_precursor_rt: Assign precursor retention times? (This relies on all precursor spectra being present and in the right order.)

        void getSpectrumMetaData(Size index, SpectrumMetaData& meta) except + nogil 
        # wrap-doc:
                #  Look up meta data of a spectrum
                #  
                #  :param index: Index of the spectrum
                #  :param meta: Meta data output

        void getSpectrumMetaData(String spectrum_ref, SpectrumMetaData& meta) except + nogil 
        # wrap-doc:
                #  Extract meta data from a spectrum
                #  
                #  :param spectrum: Spectrum input
                #  :param meta: Meta data output
                #  :param scan_regexp: Regular expression for extracting scan number from spectrum native ID
                #  :param precursor_rts: RTs of potential precursor spectra of different MS levels

        void getSpectrumMetaData(String spectrum_ref, SpectrumMetaData& meta, unsigned char flags) except + nogil 
        # wrap-doc:
                #  Extract meta data via a spectrum reference
                #  
                #  :param spectrum_ref: Spectrum reference to parse
                #  :param metadata: Meta data output
                #  :param flags: What meta data to extract

        void setSpectraDataRef(const String & spectra_data) except + nogil 

#
## wrap static methods
#
cdef extern from "<OpenMS/METADATA/SpectrumMetaDataLookup.h>" namespace "OpenMS::SpectrumMetaDataLookup":

    ###   TODO: Missing optional parameters:
    ###   boost::regex& scan_regexp 
    ###   std::map<Size, double>& precursor_rts 
    void getSpectrumMetaData(MSSpectrum spectrum,
                             SpectrumMetaData& meta) except + nogil  # wrap-attach:SpectrumMetaDataLookup

    bool addMissingRTsToPeptideIDs(libcpp_vector[PeptideIdentification], 
								  MSExperiment exp) except + nogil  # wrap-attach:SpectrumMetaDataLookup

    bool addMissingIMToPeptideIDs(libcpp_vector[PeptideIdentification], 
								  MSExperiment exp) except + nogil  # wrap-attach:SpectrumMetaDataLookup

    bool addMissingSpectrumReferences(libcpp_vector[PeptideIdentification], 
                                   String filename, bool stop_on_error, 
                                   bool override_spectra_data, 
                                   bool override_spectra_references,
                                   libcpp_vector[ProteinIdentification] proteins) except + nogil  # wrap-attach:SpectrumMetaDataLookup

cdef extern from "<OpenMS/METADATA/SpectrumMetaDataLookup.h>" namespace "OpenMS::SpectrumMetaDataLookup":

    cdef cppclass SpectrumMetaData "OpenMS::SpectrumMetaDataLookup::SpectrumMetaData":

      SpectrumMetaData() except + nogil 
      SpectrumMetaData(SpectrumMetaData &) except + nogil 

      double rt
      double precursor_rt
      double precursor_mz
      Int precursor_charge
      Size ms_level
      Int scan_number
      String native_id
