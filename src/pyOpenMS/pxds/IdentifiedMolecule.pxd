from String cimport *
from MetaData cimport *
from libcpp.set cimport set as libcpp_set
from IdentifiedCompound cimport *
from IdentifiedSequence cimport *
from EmpiricalFormula cimport *


cdef extern from "<OpenMS/METADATA/ID/IdentifiedMolecule.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass IdentifiedMolecule: #We lose anything inherited because RefVariant isn't currently supported

    IdentifiedMolecule() nogil except +

    IdentifiedMolecule(IdentifiedMolecule & other) nogil except +

    IdentifiedMolecule(IdentifiedPeptideRef ref) nogil except +
    IdentifiedMolecule(IdentifiedCompoundRef ref) nogil except +
    IdentifiedMolecule(IdentifiedOligoRef ref) nogil except +

    bool operator==(const IdentifiedMolecule & other) nogil except +
    bool operator!=(const IdentifiedMolecule & other) nogil except +
    bool operator<(const IdentifiedMolecule & other) nogil except +

    MoleculeType getMoleculeType() nogil except +

    IdentifiedPeptideRef getIdentifiedPeptideRef() nogil except +
    IdentifiedCompoundRef getIdentifiedCompoundRef() nogil except +
    IdentifiedOligoRef getIdentifiedOligoRef() nogil except +

    String toString() nogil except +

    EmpiricalFormula getFormula(Size fragment_type, int charge)
        