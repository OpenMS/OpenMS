

    def load(self, filename, peptide_ids, exp):
        """
        load(self, filename, peptide_ids, exp)
        
        Loads a map from an MSP file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the MSP file to load
        :type filename: str or bytes
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param exp: MSExperiment to store the spectra
        :type exp: MSExperiment
        
        Example::
        
            # New style (3.5+)
            peptide_ids = pyopenms.PeptideIdentificationList()
            exp = pyopenms.MSExperiment()
            pyopenms.MSPFile().load("test.msp", peptide_ids, exp)
            
            # Old style (backward compatible)
            peptide_ids = []
            exp = pyopenms.MSExperiment()
            pyopenms.MSPFile().load("test.msp", peptide_ids, exp)
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
            cpp_peptide_ids,
            deref(exp.inst.get())
        )
        
        # If we converted from Python list, copy results back
        if needs_conversion:
            peptide_ids[:] = []  # Clear the original list
            for i in range(cpp_peptide_ids.size()):
                peptide_ids.append(cpp_peptide_ids[i])
