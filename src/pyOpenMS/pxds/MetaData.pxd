from libcpp.set cimport set as libcpp_set
from String cimport *

cdef extern from "<OpenMS/METADATA/ID/MetaData.h>" namespace "OpenMS::IdentificationDataInternal":

    
    cdef cppclass IteratorWrapper[Iterator]: # not sure if this is the correct way to handle templated inheritance
        IteratorWrapper() nogil except +
        IteratorWrapper(IteratorWrapper) nogil except +
        IteratorWrapper(Iterator & it) nogil except + #why do I need this as well as the one above
        bool operator<(IteratorWrapper & other) nogil except +


    cdef enum MoleculeType:
        PROTEIN = 0,
        COMPOUND,
        RNA,
        SIZE_OF_MOLECULETYPE

    cdef enum MassType:
        MONOISOTOPIC = 0,
        AVERAGE,
        SIZE_OF_MASSTYPE

    ctypedef libcpp_set[String] InputFiles
    ctypedef IteratorWrapper[InputFiles].Iterator InputFileRef