from MzTab cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MzTabFile.h>" namespace "OpenMS":

    cdef cppclass MzTabFile:

        MzTabFile() nogil except +

        # void store(String filename, MzTab & mz_tab) nogil except +

#  -- TODO missing function in PXD:  void store(String & filename, libcpp_vector[ ProteinIdentification ] & protein_ids, libcpp_vector[ PeptideIdentification ] & peptide_ids, String in_, String document_id)
#  -- TODO missing function in PXD:  void store(String & filename, MzTab & mz_tab)
#  -- TODO missing function in PXD:  void load(String & filename, MzTab & mz_tab)
