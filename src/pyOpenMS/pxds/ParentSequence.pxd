from MetaData cimport *
from ScoredProcessingResult cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/ID/ParentSequence.h>" namespace "OpenMS::IdentificationDataInternal":
  cdef cppclass ParentSequence(ScoredProcessingResult):
    String accession

    MoleculeType molecule_type

    String sequence
    String description

    double coverage

    bool is_decoy

    ParentSequence(const String& accession, MoleculeType molecule_type, const String& sequence, const String&, double coverage, bool is_decoy , const AppliedProcessingSteps& steps_and_scores) nogil except +

    ParentSequence(const ParentSequence & other) nogil except +

    ParentSequence() nogil except +

    ParentSequence & merge(const ParentSequence & other)
        
  cdef cppclass ParentSequences:
    ParentSequences() nogil except +
    ParentSequences(ParentSequences & other) nogil except +

  cdef cppclass ParentSequenceRef:
    ParentSequenceRef() nogil except +
    ParentSequenceRef(const ParentSequenceRef & other) nogil except +
    bool operator!=(const ParentSequenceRef & other) nogil except +
    bool operator<(const ParentSequenceRef & other) nogil except +
    ParentSequence deref() nogil except +