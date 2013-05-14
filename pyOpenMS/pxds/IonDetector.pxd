from String cimport *
from Software cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/IonDetector.h>" namespace "OpenMS":

    cdef cppclass IonDetector(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        IonDetector() nogil except +
        IonDetector(IonDetector) nogil except + # wrap-ignore

        # returns the detector type
        # Type_IonDetector getType() nogil except +
        # sets the detector type
        # void setType(Type_IonDetector type_) nogil except +

        # returns the acquisition mode
        AcquisitionMode getAcquisitionMode() nogil except +
        # sets the acquisition mode
        void setAcquisitionMode(AcquisitionMode acquisition_mode) nogil except +

        # returns the resolution (in ns)
        DoubleReal getResolution() nogil except +
        # sets the resolution (in ns)
        void setResolution(DoubleReal resolution) nogil except +

        # returns the analog-to-digital converter sampling frequency (in Hz)
        DoubleReal getADCSamplingFrequency() nogil except +
        # sets the analog-to-digital converter sampling frequency (in Hz)
        void setADCSamplingFrequency(DoubleReal ADC_sampling_frequency) nogil except +

        Int getOrder() nogil except +
        # sets the order
        void setOrder(Int order) nogil except +

cdef extern from "<OpenMS/METADATA/IonDetector.h>" namespace "OpenMS::IonDetector":

        # Detector type
        # cdef enum Type:
        #   # wrap-attach:
        #   #     IonDetector
        #   TYPENULL,                                                                 #< Unknown
        #   ELECTRONMULTIPLIER,                                           #< Electron multiplier
        #   PHOTOMULTIPLIER,                                                  #< Photo multiplier
        #   FOCALPLANEARRAY,                                                  #< Focal plane array
        #   FARADAYCUP,                                                           #< Faraday cup
        #   CONVERSIONDYNODEELECTRONMULTIPLIER,           #< Conversion dynode electron multiplier
        #   CONVERSIONDYNODEPHOTOMULTIPLIER,                  #< Conversion dynode photo multiplier
        #   MULTICOLLECTOR,                                                   #< Multi-collector
        #   CHANNELELECTRONMULTIPLIER,                            #< Channel electron multiplier
        #   CHANNELTRON,                                                          #< channeltron
        #   DALYDETECTOR,                                                         #< daly detector
        #   MICROCHANNELPLATEDETECTOR,                            #< microchannel plate detector
        #   ARRAYDETECTOR,                                                    #< array detector
        #   CONVERSIONDYNODE,                                                 #< conversion dynode
        #   DYNODE,                                                                   #< dynode
        #   FOCALPLANECOLLECTOR,                                          #< focal plane collector
        #   IONTOPHOTONDETECTOR,                                          #< ion-to-photon detector
        #   POINTCOLLECTOR,                                                   #< point collector
        #   POSTACCELERATIONDETECTOR,                                 #< postacceleration detector
        #   PHOTODIODEARRAYDETECTOR,                                  #< photodiode array detector
        #   INDUCTIVEDETECTOR,                                            #< inductive detector
        #   ELECTRONMULTIPLIERTUBE,                                   #< electron multiplier tube
        #   SIZE_OF_TYPE

        # Acquisition mode
        cdef enum AcquisitionMode:
          # wrap-attach:
          #     IonDetector
          ACQMODENULL,                          #< Unknown
          PULSECOUNTING,                    #< Pulse counting
          ADC,                                          #< Analog-digital converter
          TDC,                                          #< Time-digital converter
          TRANSIENTRECORDER,            #< Transient recorder
          SIZE_OF_ACQUISITIONMODE
