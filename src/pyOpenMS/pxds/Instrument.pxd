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
        #    MetaInfoInterface

        Instrument() nogil except + # wrap-doc:Description of a MS instrument
        Instrument(Instrument &) nogil except +

        # returns the name of the instrument
        String getName() nogil except + # wrap-doc:Returns the name of the instrument
        # sets the name of the instrument
        void setName(String name) nogil except + # wrap-doc:Sets the name of the instrument

        # returns the instrument vendor
        String getVendor() nogil except + # wrap-doc:Returns the instrument vendor
        # sets the instrument vendor
        void setVendor(String vendor) nogil except + # wrap-doc:Sets the instrument vendor

        # returns the instrument model
        String getModel() nogil except + # wrap-doc:Returns the instrument model
        # sets the instrument model
        void setModel(String model) nogil except + # wrap-doc:Sets the instrument model

        # returns a description of customizations
        String getCustomizations() nogil except + # wrap-doc:Returns a description of customizations
        # sets the a description of customizations
        void setCustomizations(String customizations) nogil except + # wrap-doc:Sets the a description of customizations

        # returns a mutable reference to the ion source list
        libcpp_vector[IonSource] getIonSources() nogil except + # wrap-doc:Returns the ion source list
        # sets the ion source list
        void setIonSources(libcpp_vector[IonSource] ion_sources) nogil except + # wrap-doc:Sets the ion source list

        # returns a mutable reference to the mass analyzer list
        libcpp_vector[MassAnalyzer] getMassAnalyzers() nogil except + # wrap-doc:Returns the mass analyzer list
        # sets the mass analyzer list
        void setMassAnalyzers(libcpp_vector[MassAnalyzer] mass_analyzers) nogil except + # wrap-doc:Sets the mass analyzer list

        # returns a mutable reference to the ion detector list
        libcpp_vector[IonDetector] getIonDetectors() nogil except + # wrap-doc:Returns the ion detector list
        # sets the ion detector list
        void setIonDetectors(libcpp_vector[IonDetector] ion_detectors) nogil except + # wrap-doc:Sets the ion detector list

        # returns a mutable reference to the instrument software
        Software getSoftware() nogil except + # wrap-doc:Returns the instrument software
        # sets the instrument software
        void setSoftware(Software software) nogil except + # wrap-doc:Sets the instrument software

        IonOpticsType getIonOptics() nogil except + # wrap-doc:Returns the ion optics type
        void setIonOptics(IonOpticsType ion_optics) nogil except + # wrap-doc:Sets the ion optics type

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

