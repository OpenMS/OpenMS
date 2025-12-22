

    def load(self, filename, protein_ids, peptide_ids, *args):
        """
        load(self, filename, protein_ids, peptide_ids, experiment_name=None, lookup=None)
        
        Loads the identifications from a PepXML file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the PepXML file to load
        :type filename: str or bytes
        :param protein_ids: List to store protein identifications (modified in place)
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param experiment_name: Optional experiment name
        :type experiment_name: str or bytes, optional
        :param lookup: Optional spectrum metadata lookup structure
        :type lookup: SpectrumMetaDataLookup, optional
        
        Example::
        
            # New style (3.5+)
            protein_ids = []
            peptide_ids = pyopenms.PeptideIdentificationList()
            pyopenms.PepXMLFile().load("test.pep.xml", protein_ids, peptide_ids)
            
            # Old style (backward compatible)
            protein_ids = []
            peptide_ids = []
            pyopenms.PepXMLFile().load("test.pep.xml", protein_ids, peptide_ids)
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
        
        # Call the appropriate C++ method based on number of arguments
        if len(args) == 0:
            # load(filename, protein_ids, peptide_ids)
            self.inst.get().load(
                _String(<char*>_bytes(filename)),
                protein_ids,
                cpp_peptide_ids
            )
        elif len(args) == 1:
            # load(filename, protein_ids, peptide_ids, experiment_name)
            experiment_name = args[0]
            self.inst.get().load(
                _String(<char*>_bytes(filename)),
                protein_ids,
                cpp_peptide_ids,
                _String(<char*>_bytes(experiment_name))
            )
        elif len(args) == 2:
            # load(filename, protein_ids, peptide_ids, experiment_name, lookup)
            experiment_name = args[0]
            lookup = args[1]
            self.inst.get().load(
                _String(<char*>_bytes(filename)),
                protein_ids,
                cpp_peptide_ids,
                _String(<char*>_bytes(experiment_name)),
                deref(lookup.inst.get())
            )
        else:
            raise ValueError("Invalid number of arguments")
        
        # If we converted from Python list, copy results back
        if needs_conversion:
            peptide_ids[:] = []  # Clear the original list
            for i in range(cpp_peptide_ids.size()):
                peptide_ids.append(cpp_peptide_ids[i])
    
    def store(self, filename, protein_ids, peptide_ids, *args):
        """
        store(self, filename, protein_ids, peptide_ids, mz_file="", mz_name="", 
              peptideprophet_analyzed=False, rt_tolerance=0.0)
        
        Stores the identifications to a PepXML file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the PepXML file to store
        :type filename: str or bytes
        :param protein_ids: List of protein identifications to store
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container of peptide identifications to store.
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param mz_file: Optional mz file name
        :type mz_file: str or bytes, optional
        :param mz_name: Optional mz name
        :type mz_name: str or bytes, optional
        :param peptideprophet_analyzed: Whether data was analyzed with PeptideProphet
        :type peptideprophet_analyzed: bool, optional
        :param rt_tolerance: RT tolerance
        :type rt_tolerance: float, optional
        
        Example::
        
            # New style (3.5+)
            protein_ids = [...]
            peptide_ids = pyopenms.PeptideIdentificationList()
            peptide_ids.extend([...])
            pyopenms.PepXMLFile().store("test.pep.xml", protein_ids, peptide_ids)
            
            # Old style (backward compatible)
            protein_ids = [...]
            peptide_ids = [...]
            pyopenms.PepXMLFile().store("test.pep.xml", protein_ids, peptide_ids)
        """
        cdef PeptideIdentificationList cpp_peptide_ids
        
        # Check if peptide_ids is a Python list (old API) or PeptideIdentificationList (new API)
        if isinstance(peptide_ids, list):
            # Old API: Convert Python list to PeptideIdentificationList
            cpp_peptide_ids = PeptideIdentificationList()
            for pep_id in peptide_ids:
                cpp_peptide_ids.push_back(pep_id)
        else:
            # New API: Use PeptideIdentificationList directly
            cpp_peptide_ids = peptide_ids
        
        # Call the appropriate C++ method based on number of arguments
        if len(args) == 0:
            # store(filename, protein_ids, peptide_ids)
            self.inst.get().store(
                _String(<char*>_bytes(filename)),
                protein_ids,
                cpp_peptide_ids
            )
        elif len(args) == 4:
            # store(filename, protein_ids, peptide_ids, mz_file, mz_name, peptideprophet_analyzed, rt_tolerance)
            mz_file = args[0]
            mz_name = args[1]
            peptideprophet_analyzed = args[2]
            rt_tolerance = args[3]
            self.inst.get().store(
                _String(<char*>_bytes(filename)),
                protein_ids,
                cpp_peptide_ids,
                _String(<char*>_bytes(mz_file)),
                _String(<char*>_bytes(mz_name)),
                peptideprophet_analyzed,
                rt_tolerance
            )
        else:
            raise ValueError("Invalid number of arguments")
