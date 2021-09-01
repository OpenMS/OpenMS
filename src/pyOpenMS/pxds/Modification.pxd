from Types cimport *
from SampleTreatment cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/Modification.h>" namespace "OpenMS":
    
    cdef cppclass Modification(SampleTreatment) :
        # wrap-inherits:
        #  SampleTreatment
        Modification() nogil except +
        Modification(Modification &) nogil except +

        # bool operator==(SampleTreatment & rhs) nogil except +

        # SampleTreatment * clone() nogil except +

        String getReagentName() nogil except + # wrap-doc:Returns the name of the reagent that was used (default "")
        void setReagentName(const String & reagent_name) nogil except + # wrap-doc:Sets the name of the reagent that was used

        double getMass() nogil except + # wrap-doc:Returns the mass change (default 0.0)
        void setMass(double mass) nogil except + # wrap-doc:Sets the mass change

        Modification_SpecificityType getSpecificityType() nogil except + # wrap-doc:Returns the specificity of the reagent (default AA)
        void setSpecificityType(Modification_SpecificityType & specificity_type) nogil except + # wrap-doc:Sets the specificity of the reagent

        String getAffectedAminoAcids() nogil except + # wrap-doc:Returns a string containing the one letter code of the amino acids that are affected by the reagent (default "")
        void setAffectedAminoAcids(const String & affected_amino_acids) nogil except + # wrap-doc:Returns a string containing the one letter code of the amino acids that are affected by the reagent. Do not separate them by space, tab or comma!

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
