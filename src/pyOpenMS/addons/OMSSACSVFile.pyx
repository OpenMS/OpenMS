

    def load(self, filename, protein_identification, peptide_ids):
        """
        load(self, filename, protein_identification, peptide_ids)
        
        Loads data from an OMSSA CSV file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the OMSSA CSV file to load
        :type filename: str or bytes
        :param protein_identification: Protein identification object to store results (modified in place)
        :type protein_identification: ProteinIdentification
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        
        Example::
        
            # New style (3.5+)
            protein_id = pyopenms.ProteinIdentification()
            peptide_ids = pyopenms.PeptideIdentificationList()
            pyopenms.OMSSACSVFile().load("test.csv", protein_id, peptide_ids)
            
            # Old style (backward compatible)
            protein_id = pyopenms.ProteinIdentification()
            peptide_ids = []
            pyopenms.OMSSACSVFile().load("test.csv", protein_id, peptide_ids)
        """
        cdef PeptideIdentificationList cpp_peptide_ids
        cdef bool needs_conversion = False
        
        # Check if peptide_ids is a Python list (old API) or PeptideIdentificationList (new API)
        if isinstance(peptide_ids, list):
            # Old API: Convert Python list to PeptideIdentificationList
            needs_conversion = True
            cpp_peptide_ids = PeptideIdentificationList()
            # Copy existing items from Python list to PeptideIdentificationList
            for pep_id in peptide_ids:
                cpp_peptide_ids.push_back(pep_id)
        else:
            # New API: Use PeptideIdentificationList directly
            cpp_peptide_ids = peptide_ids
        
        # Call the C++ load method
        self.inst.get().load(
            _String(<char*>_bytes(filename)),
            deref(protein_identification.inst.get()),
            cpp_peptide_ids
        )
        
        # If we converted from Python list, copy results back
        if needs_conversion:
            peptide_ids[:] = []  # Clear the original list
            for i in range(cpp_peptide_ids.size()):
                peptide_ids.append(cpp_peptide_ids[i])
