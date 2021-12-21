from Types cimport *
from String cimport String

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":
    
    cdef cppclass CsiAdapterHit "OpenMS::CsiFingerIdMzTabWriter::CsiAdapterHit":

        CsiAdapterHit() nogil except + 
        CsiAdapterHit(CsiAdapterHit &) nogil except + # compiler

        String inchikey2D
        String inchi
        unsigned int rank
        unsigned int formula_rank
        String adduct
        String molecular_formula
        double score
        String name
        String smiles
        String xlogp
        String dbflags
        libcpp_vector[ String ] pubchemids
        libcpp_vector[ String ] links
