from MetaInfoInterface cimport *
from IonSource cimport *
from ScanWindow cimport *

cdef extern from "<OpenMS/METADATA/InstrumentSettings.h>" namespace "OpenMS":

    cdef cppclass InstrumentSettings(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        InstrumentSettings()     nogil except +
        InstrumentSettings(InstrumentSettings)     nogil except +

        Polarity getPolarity()     nogil except +
        void setPolarity(Polarity)  nogil except +

        ScanMode getScanMode() nogil except +
        void setScanMode(ScanMode scan_mode) nogil except +
        bool getZoomScan() nogil except +
        void setZoomScan(bool zoom_scan) nogil except +
        libcpp_vector[ ScanWindow ]  getScanWindows() nogil except +
        void setScanWindows(libcpp_vector[ ScanWindow ] scan_windows) nogil except +

        # declare again: cython complains for overloaded methods in base
        # classes
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +


cdef extern from "<OpenMS/METADATA/InstrumentSettings.h>" namespace "OpenMS::InstrumentSettings":

    # scan mode
    cdef enum ScanMode:
      UNKNOWN,                      #< Unknown scan method
      MASSSPECTRUM,       #< general spectrum type
      MS1SPECTRUM,              #< full scan mass spectrum, is a "mass spectrum" @n Synonyms: 'full spectrum', 'Q1 spectrum', 'Q3 spectrum', 'Single-Stage Mass Spectrometry'
      MSNSPECTRUM,        #< MS2+ mass spectrum, is a "mass spectrum"
      SIM,                              #< Selected ion monitoring scan @n Synonyms: 'Multiple ion monitoring scan', 'SIM scan', 'MIM scan'
      SRM,                              #< Selected reaction monitoring scan @n Synonyms: 'Multiple reaction monitoring scan', 'SRM scan', 'MRM scan'
      CRM,                              #< Consecutive reaction monitoring scan @n Synonyms: 'CRM scan'
      CNG,                              #< Constant neutral gain scan @n Synonyms: 'CNG scan'
      CNL,                              #< Constant neutral loss scan @n Synonyms: 'CNG scan'
      PRECURSOR,                #< Precursor ion scan
      EMC,                              #< Enhanced multiply charged scan
      TDF,                              #< Time-delayed fragmentation scan
      EMR,                              #< Electromagnetic radiation scan @n Synonyms: 'EMR spectrum'
      EMISSION,                     #< Emission scan
      ABSORPTION,               #< Absorption scan
      SIZE_OF_SCANMODE

