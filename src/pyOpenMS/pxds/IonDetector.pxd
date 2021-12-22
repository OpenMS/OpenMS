from String cimport *
from Software cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/IonDetector.h>" namespace "OpenMS":

    cdef cppclass IonDetector(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface
        
        IonDetector() nogil except + # wrap-doc:Description of a ion detector (part of a MS Instrument)
        IonDetector(IonDetector &) nogil except + 

        Type_IonDetector getType() nogil except + # wrap-doc:Returns the detector type
        void setType(Type_IonDetector type_) nogil except + # wrap-doc:Sets the detector type

        AcquisitionMode getAcquisitionMode() nogil except + # wrap-doc:Returns the acquisition mode
        void setAcquisitionMode(AcquisitionMode acquisition_mode) nogil except + # wrap-doc:Sets the acquisition mode

        double getResolution() nogil except + # wrap-doc:Returns the resolution (in ns)
        void setResolution(double resolution) nogil except + # wrap-doc:Sets the resolution (in ns)

        double getADCSamplingFrequency() nogil except + # wrap-doc:Returns the analog-to-digital converter sampling frequency (in Hz)
        void setADCSamplingFrequency(double ADC_sampling_frequency) nogil except + # wrap-doc:Sets the analog-to-digital converter sampling frequency (in Hz)

        Int getOrder() nogil except + # wrap-doc:Returns the order
        void setOrder(Int order) nogil except + # wrap-doc:Sets the order

cdef extern from "<OpenMS/METADATA/IonDetector.h>" namespace "OpenMS::IonDetector":

        # Detector type
        cdef enum Type_IonDetector "OpenMS::IonDetector::Type":
          # wrap-attach:
          #     IonDetector
          TYPENULL,                                                                 #< Unknown
          ELECTRONMULTIPLIER,                                           #< Electron multiplier
          PHOTOMULTIPLIER,                                                  #< Photo multiplier
          FOCALPLANEARRAY,                                                  #< Focal plane array
          FARADAYCUP,                                                           #< Faraday cup
          CONVERSIONDYNODEELECTRONMULTIPLIER,           #< Conversion dynode electron multiplier
          CONVERSIONDYNODEPHOTOMULTIPLIER,                  #< Conversion dynode photo multiplier
          MULTICOLLECTOR,                                                   #< Multi-collector
          CHANNELELECTRONMULTIPLIER,                            #< Channel electron multiplier
          CHANNELTRON,                                                          #< channeltron
          DALYDETECTOR,                                                         #< daly detector
          MICROCHANNELPLATEDETECTOR,                            #< microchannel plate detector
          ARRAYDETECTOR,                                                    #< array detector
          CONVERSIONDYNODE,                                                 #< conversion dynode
          DYNODE,                                                                   #< dynode
          FOCALPLANECOLLECTOR,                                          #< focal plane collector
          IONTOPHOTONDETECTOR,                                          #< ion-to-photon detector
          POINTCOLLECTOR,                                                   #< point collector
          POSTACCELERATIONDETECTOR,                                 #< postacceleration detector
          PHOTODIODEARRAYDETECTOR,                                  #< photodiode array detector
          INDUCTIVEDETECTOR,                                            #< inductive detector
          ELECTRONMULTIPLIERTUBE,                                   #< electron multiplier tube
          SIZE_OF_TYPE

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
