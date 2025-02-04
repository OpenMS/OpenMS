from String cimport *
from IonSource cimport *
from MassAnalyzer cimport *
from Software cimport *
from IonDetector cimport *
from MetaInfoInterface cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/Instrument.h>" namespace "OpenMS":

    cdef cppclass Instrument(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface
        # wrap-doc:
        #  Description of a MS instrument

        Instrument() except + nogil  # wrap-doc:Description of a MS instrument
        Instrument(Instrument &) except + nogil 

        String getName() except + nogil  # wrap-doc:Returns the name of the instrument
        void setName(String name) except + nogil  # wrap-doc:Sets the name of the instrument

        String getVendor() except + nogil  # wrap-doc:Returns the instrument vendor
        void setVendor(String vendor) except + nogil  # wrap-doc:Sets the instrument vendor

        String getModel() except + nogil  # wrap-doc:Returns the instrument model
        void setModel(String model) except + nogil  # wrap-doc:Sets the instrument model

        String getCustomizations() except + nogil  # wrap-doc:Returns a description of customizations
        void setCustomizations(String customizations) except + nogil  # wrap-doc:Sets the a description of customizations

        libcpp_vector[IonSource] getIonSources() except + nogil  # wrap-doc:Returns the ion source list
        void setIonSources(libcpp_vector[IonSource] ion_sources) except + nogil  # wrap-doc:Sets the ion source list

        libcpp_vector[MassAnalyzer] getMassAnalyzers() except + nogil  # wrap-doc:Returns the mass analyzer list
        void setMassAnalyzers(libcpp_vector[MassAnalyzer] mass_analyzers) except + nogil  # wrap-doc:Sets the mass analyzer list

        libcpp_vector[IonDetector] getIonDetectors() except + nogil  # wrap-doc:Returns the ion detector list
        void setIonDetectors(libcpp_vector[IonDetector] ion_detectors) except + nogil  # wrap-doc:Sets the ion detector list

        Software getSoftware() except + nogil  # wrap-doc:Returns the instrument software
        void setSoftware(Software software) except + nogil  # wrap-doc:Sets the instrument software

        IonOpticsType getIonOptics() except + nogil  # wrap-doc:Returns the ion optics type
        void setIonOptics(IonOpticsType ion_optics) except + nogil  # wrap-doc:Sets the ion optics type

cdef extern from "<OpenMS/METADATA/Instrument.h>" namespace "OpenMS::Instrument":

    cdef enum IonOpticsType:
    
        UNKNOWN,                          #< unknown
        MAGNETIC_DEFLECTION,              #< magnetic deflection
        DELAYED_EXTRACTION,               #< delayed extraction
        COLLISION_QUADRUPOLE,             #< collision quadrupole
        SELECTED_ION_FLOW_TUBE,           #< selected ion flow tube
        TIME_LAG_FOCUSING,                #< time lag focusing
        REFLECTRON,                       #< reflectron
        EINZEL_LENS,                      #< einzel lens
        FIRST_STABILITY_REGION,           #< first stability region
        FRINGING_FIELD,                   #< fringing field
        KINETIC_ENERGY_ANALYZER,          #< kinetic energy analyzer
        STATIC_FIELD,                     #< static field
        SIZE_OF_IONOPTICSTYPE
