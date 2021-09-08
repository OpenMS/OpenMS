from libcpp.vector cimport vector as libcpp_vector
from String cimport *

from ProteinIdentification cimport *
from PeptideIdentification cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

from libcpp cimport bool

cdef extern from "<OpenMS/FORMAT/ProtXMLFile.h>" namespace "OpenMS":

    cdef cppclass ProtXMLFile:
        # wrap-doc:
            #   Used to load (storing not supported, yet) ProtXML files
            #   -----
            #   This class is used to load (storing not supported, yet) documents that implement
            #   the schema of ProtXML files
            #
            #   A documented schema for this format comes with the TPP and can also be
            #   found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS
            #   -----
            #   OpenMS can only read parts of the protein_summary subtree to extract
            #   protein-peptide associations. All other parts are silently ignored
            #   -----
            #   For protein groups, only the "group leader" (which is annotated with a
            #   probability and coverage) receives these attributes. All indistinguishable
            #   proteins of the same group only have an accession and score of -1

        ProtXMLFile() nogil except +
        # copy constructor of 'ProtXMLFile' is implicitly deleted because base class 'Internal::XMLHandler' has a deleted copy constructor protected Internal::XMLHandler
        ProtXMLFile(ProtXMLFile &) nogil except + # wrap-ignore

        void load(String filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids) nogil except +
            # wrap-doc:
            #   Loads the identifications of an ProtXML file without identifier
            #   -----
            #   The information is read in and the information is stored in the
            #   corresponding variables
            #   -----
            #   :raises:
            #     Exception: FileNotFound is thrown if the file could not be found
            #   :raises:
            #     Exception: ParseError is thrown if an error occurs during parsing

        # Not implemented
        void store(String filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids, String document_id) nogil except +

