from String cimport *
from Software cimport *
from MetaInfoInterface cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/IonSource.h>" namespace "OpenMS":

    cdef cppclass IonSource(MetaInfoInterface):
        # wrap-inherits:
        #    MetaInfoInterface

        IonSource() nogil except +
        IonSource(IonSource) nogil except + # wrap-ignore

        # returns the ionization mode
        Polarity getPolarity() nogil except +
        # sets the ionization mode
        void setPolarity(Polarity polarity) nogil except +

        # returns the inlet type
        InletType getInletType() nogil except +
        # sets the  inlet type
        void setInletType(InletType inlet_type) nogil except +

        # returns the ionization method
        IonizationMethod getIonizationMethod() nogil except +
        # sets the ionization method
        void setIonizationMethod(IonizationMethod ionization_type) nogil except +

        #     @brief returns the position of this part in the whole Instrument.
        #     Order can be ignored, as long the instrument has this default setup:
        #     - one ion source
        #     - one or many mass analyzers
        #     - one ion detector
        #     For more complex instruments, the order should be defined.
        Int getOrder() nogil except +
        # sets the order
        void setOrder(Int order) nogil except +

cdef extern from "<OpenMS/METADATA/IonSource.h>" namespace "OpenMS::IonSource":

    cdef enum Polarity:
        # wrap-attach:
        #     IonSource
        POLNULL, POSITIVE, NEGATIVE, SIZE_OF_POLARITY

    cdef enum InletType:
        # wrap-attach:
        #     IonSource
        INLETNULL,                                                        #]Unknown
        DIRECT,                                                               #]Direct
        BATCH,                                                                #]Batch (e.g. in MALDI)
        CHROMATOGRAPHY,                                               #]Chromatography (liquid)
        PARTICLEBEAM,                                                     #]Particle beam
        MEMBRANESEPARATOR,                                        #]Membrane sparator
        OPENSPLIT,                                                        #]Open split
        JETSEPARATOR,                                                     #]Jet separator
        SEPTUM,                                                               #]Septum
        RESERVOIR,                                                        #]Reservoir
        MOVINGBELT,                                                       #]Moving belt
        MOVINGWIRE,                                                       #]Moving wire
        FLOWINJECTIONANALYSIS,                                #]Flow injection analysis
        ELECTROSPRAYINLET,                                        #]Electro spray
        THERMOSPRAYINLET,                                             #]Thermo spray
        INFUSION,                                                             #]Infusion
        CONTINUOUSFLOWFASTATOMBOMBARDMENT,        #]Continuous flow fast atom bombardment
        INDUCTIVELYCOUPLEDPLASMA,                             #]Inductively coupled plasma
        MEMBRANE,                                                             #]Membrane inlet
        NANOSPRAY,                                                        #]Nanospray inlet
        SIZE_OF_INLETTYPE


    # ionization method
    cdef enum IonizationMethod:
        # wrap-attach:
        #     IonSource
        IONMETHODNULL,        #]Unknown
        ESI,                              #]electrospray ionisation
        EI,                               #]electron ionization
        CI,                               #]chemical ionisation
        FAB,                              #]fast atom bombardment
        TSP,                              #]thermospray
        LD,                               #]laser desorption
        FD,                               #]field desorption
        FI,                               #]flame ionization
        PD,                               #]plasma desorption
        SI,                               #]secondary ion MS
        TI,                               #]thermal ionization
        API,                              #]atmospheric pressure ionisation
        ISI,                              #<
        CID,                              #]collsion induced decomposition
        CAD,                              #]collsion activated decomposition
        HN,                               #<
        APCI,                             #]atmospheric pressure chemical ionization
        APPI,                             #]atmospheric pressure photo ionization
        ICP,                              #]inductively coupled plasma
        NESI,                             #]Nano electrospray ionization
        MESI,                             #]Micro electrospray ionization
        SELDI,                        #]Surface enhanced laser desorption ionization
        SEND,                             #]Surface enhanced neat desorption
        FIB,                              #]Fast ion bombardment
        MALDI,                        #]Matrix-assisted laser desorption ionization
        MPI,                              #]Multiphoton ionization
        DI,                               #]desorption ionization
        FA,                               #]flowing afterglow
        FII,                              #]field ionization
        GD_MS,                        #]glow discharge ionization
        NICI,                             #]negative ion chemical ionization
        NRMS,                             #]neutralization reionization mass spectrometry
        PI,                               #]photoionization
        PYMS,                             #]pyrolysis mass spectrometry
        REMPI,                        #]resonance enhanced multiphoton ionization
        AI,                               #]adiabatic ionization
        ASI,                              #]associative ionization
        AD,                               #]autodetachment
        AUI,                              #]autoionization
        CEI,                              #]charge exchange ionization
        CHEMI,                        #]chemi-ionization
        DISSI,                        #]dissociative ionization
        LSI,                              #]liquid secondary ionization
        PEI,                              #]penning ionization
        SOI,                              #]soft ionization
        SPI,                              #]spark ionization
        SUI,                              #]surface ionization
        VI,                               #]vertical ionization
        AP_MALDI,                     #]atmospheric pressure matrix-assisted laser desorption ionization
        SILI,                             #]desorption/ionization on silicon
        SALDI,                        #]surface-assisted laser desorption ionization
        SIZE_OF_IONIZATIONMETHOD
