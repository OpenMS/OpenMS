from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from IdentificationData cimport *
from FASTAFile cimport *
from MzTab cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/METADATA/ID/IdentificationDataConverter.h>" namespace "OpenMS":
    
    cdef cppclass IdentificationDataConverter "OpenMS::IdentificationDataConverter":
        IdentificationDataConverter() nogil except +
        IdentificationDataConverter(IdentificationDataConverter) nogil except + #wrap-ignore
        void importIDs(IdentificationData & id_data, const libcpp_vector[ ProteinIdentification ] & proteins, const libcpp_vector[ PeptideIdentification ] & peptides) nogil except +
        void exportIDs(const IdentificationData & id_data, libcpp_vector[ ProteinIdentification ] & proteins, libcpp_vector[ PeptideIdentification ] & peptides, bool export_oligonucleotides) nogil except +
        MzTab exportMzTab(const IdentificationData & id_data) nogil except +
        #void importSequences(IdentificationData & id_data, const libcpp_vector[ FASTAFile::FASTAEntry ] & fasta, IdentificationData::MoleculeType type_, const String & decoy_pattern) nogil except +

