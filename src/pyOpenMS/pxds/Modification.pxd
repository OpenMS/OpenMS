from Types cimport *
from SampleTreatment cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/Modification.h>" namespace "OpenMS":
    
    cdef cppclass Modification(SampleTreatment) :
        # wrap-inherits:
        #  SampleTreatment
        Modification() nogil except +
        Modification(Modification) nogil except +

        # bool operator==(SampleTreatment & rhs) nogil except +

        # SampleTreatment * clone() nogil except +

        String getReagentName() nogil except +
        void setReagentName(const String & reagent_name) nogil except +

        double getMass() nogil except +
        void setMass(double mass) nogil except +

        Modification_SpecificityType getSpecificityType() nogil except +
        void setSpecificityType(Modification_SpecificityType & specificity_type) nogil except +

        String getAffectedAminoAcids() nogil except +
        void setAffectedAminoAcids(const String & affected_amino_acids) nogil except +

cdef extern from "<OpenMS/METADATA/Modification.h>" namespace "OpenMS::Modification":
    cdef enum Modification_SpecificityType "OpenMS::Modification::SpecificityType":
        #wrap-attach:
        #    Modification
        AA
        AA_AT_CTERM
        AA_AT_NTERM
        CTERM
        NTERM
        SIZE_OF_SPECIFICITYTYPE

