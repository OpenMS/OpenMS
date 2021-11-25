from MetaInfoInterface cimport *
from IonSource cimport *
from ScanWindow cimport *

cdef extern from "<OpenMS/METADATA/InstrumentSettings.h>" namespace "OpenMS":

    cdef cppclass InstrumentSettings(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        InstrumentSettings() nogil except + # wrap-doc:Description of the settings a MS Instrument was run with
        InstrumentSettings(InstrumentSettings &) nogil except +

        Polarity getPolarity()     nogil except + # wrap-doc:Returns the polarity
        void setPolarity(Polarity)  nogil except + # wrap-doc:Sets the polarity

        ScanMode getScanMode() nogil except + # wrap-doc:Returns the scan mode
        void setScanMode(ScanMode scan_mode) nogil except + # wrap-doc:Sets the scan mode
        bool getZoomScan() nogil except + # wrap-doc:Returns if this scan is a zoom (enhanced resolution) scan
        void setZoomScan(bool zoom_scan) nogil except + # wrap-doc:Sets if this scan is a zoom (enhanced resolution) scan
        libcpp_vector[ ScanWindow ]  getScanWindows() nogil except + # wrap-doc:Returns the m/z scan windows
        void setScanWindows(libcpp_vector[ ScanWindow ] scan_windows) nogil except + # wrap-doc:Sets the m/z scan windows

cdef extern from "<OpenMS/METADATA/InstrumentSettings.h>" namespace "OpenMS::InstrumentSettings":

    # scan mode
    cdef enum ScanMode:
      UNKNOWN,                #< Unknown scan method
      MASSSPECTRUM,           #< general spectrum type
      MS1SPECTRUM,            #< full scan mass spectrum, is a "mass spectrum" @n Synonyms: 'full spectrum', 'Q1 spectrum', 'Q3 spectrum', 'Single-Stage Mass Spectrometry'
      MSNSPECTRUM,            #< MS2+ mass spectrum, is a "mass spectrum"
      SIM,                    #< Selected ion monitoring scan @n Synonyms: 'Multiple ion monitoring scan', 'SIM scan', 'MIM scan'
      SRM,                    #< Selected reaction monitoring scan @n Synonyms: 'Multiple reaction monitoring scan', 'SRM scan', 'MRM scan'
      CRM,                    #< Consecutive reaction monitoring scan @n Synonyms: 'CRM scan'
      CNG,                    #< Constant neutral gain scan @n Synonyms: 'CNG scan'
      CNL,                    #< Constant neutral loss scan @n Synonyms: 'CNG scan'
      PRECURSOR,              #< Precursor ion scan
      EMC,                    #< Enhanced multiply charged scan
      TDF,                    #< Time-delayed fragmentation scan
      EMR,                    #< Electromagnetic radiation scan @n Synonyms: 'EMR spectrum'
      EMISSION,               #< Emission scan
      ABSORPTION,             #< Absorption scan
      SIZE_OF_SCANMODE
