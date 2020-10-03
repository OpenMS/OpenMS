from Types cimport *
from String cimport *
from MzTab cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS":
    
    cdef cppclass CsiFingerIdMzTabWriter "OpenMS::CsiFingerIdMzTabWriter":
        CsiFingerIdMzTabWriter() nogil except +  # wrap-ignore
        CsiFingerIdMzTabWriter(CsiFingerIdMzTabWriter) nogil except + #wrap-ignore

#
# wrap static method:
#
cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":

   void read(libcpp_vector[String] sirius_output_paths, String original_input_mzml, Size top_n_hits, MzTab & result) nogil except + # wrap-attach:CsiFingerIdMzTabWriter

cdef extern from "<OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>" namespace "OpenMS::CsiFingerIdMzTabWriter":

    cdef cppclass CsiAdapterIdentification "OpenMS::CsiFingerIdMzTabWriter::CsiAdapterIdentification":
        CsiAdapterIdentification() nogil except +
        CsiAdapterIdentification(CsiAdapterIdentification) nogil except + # wrap-ignore

        double mz
        double rt
        String native_id
        int scan_index
        int scan_number
        String feature_id
        libcpp_vector[CsiAdapterHit] hits

    cdef cppclass CsiAdapterHit "OpenMS::CsiFingerIdMzTabWriter::CsiAdapterHit":
        CsiAdapterHit() nogil except +
        CsiAdapterHit(CsiAdapterHit) nogil except + # wrap-ignore

        String inchikey2D
        String inchi
        int rank
        String molecular_formula
        double score
        String name
        String smiles
        libcpp_vector[String] pubchemids
        libcpp_vector[String] links
