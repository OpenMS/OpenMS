from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from StringList cimport *
from AASequence cimport *
from NASequence cimport *

cdef extern from "<OpenMS/METADATA/ID/MetaData.h>" namespace "OpenMS::IdentificationDataInternal::MoleculeType":

    ctypedef enum MoleculeType "OpenMS::IdentificationDataInternal::MoleculeType":
            PROTEIN,
            COMPOUND,
            RNA,
            SIZE_OF_MOLECULETYPE

cdef extern from "<OpenMS/METADATA/ID/ParentMolecule.h>" namespace "OpenMS::IdentificationDataInternal":
    
    cdef cppclass ParentMolecule "OpenMS::IdentificationDataInternal::ParentMolecule":
        ParentMolecule() nogil except + # wrap-ignore
        ParentMolecule(const String& accession) nogil except +
        ParentMolecule(
            const String& accession,
            MoleculeType molecule_type,
            const String& sequence, const String& description,
            double coverage, bool is_decoy) nogil except + 

        ParentMolecule(const ParentMolecule&) nogil except +

        # members
        String accession
        MoleculeType molecule_type
        String sequence
        String description
        double coverage
        bool is_decoy;

    cdef cppclass ParentMoleculeRef "OpenMS::IdentificationDataInternal::ParentMoleculeRef":
        # ParentMoleculeRef() nogil except + # wrap-ignore
        ParentMoleculeRef(ParentMoleculeRef&) nogil except + # compiler

        ParentMolecule getParentMolecule() nogil except +

        bool operator==(ParentMoleculeRef) nogil except +
        # bool operator<ParentMoleculeRef) nogil except +
        # bool operator!=(ParentMoleculeRef) nogil except +

