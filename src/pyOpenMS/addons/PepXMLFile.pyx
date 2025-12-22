

    def load(self, filename, protein_ids, peptide_ids, experiment_name=None, lookup=None):
        """
        load(self, filename, protein_ids, peptide_ids, experiment_name=None, lookup=None)
        
        Loads the identifications from a PepXML file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the PepXML file to load
        :type filename: str or bytes or String
        :param protein_ids: List to store protein identifications (modified in place)
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param experiment_name: Optional experiment name
        :type experiment_name: str or bytes or String, optional
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
        # Import at function level to avoid circular dependencies
        from . import PeptideIdentificationList as _PeptideIdentificationList
        
        # Check if peptide_ids is a Python list (old API) or PeptideIdentificationList (new API)
        if isinstance(peptide_ids, list):
            # Old API: Convert Python list to PeptideIdentificationList
            temp_peptide_ids = _PeptideIdentificationList()
            temp_peptide_ids.extend(peptide_ids)
            
            # Call the appropriate C++ method based on arguments
            if experiment_name is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(temp_peptide_ids.inst.get())
                )
            elif lookup is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(temp_peptide_ids.inst.get()),
                    deref(convString(experiment_name).get())
                )
            else:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(temp_peptide_ids.inst.get()),
                    deref(convString(experiment_name).get()),
                    deref(lookup.inst.get())
                )
            
            # Copy results back to Python list
            peptide_ids[:] = list(temp_peptide_ids)
        else:
            # New API: Use PeptideIdentificationList directly
            if experiment_name is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(peptide_ids.inst.get())
                )
            elif lookup is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(peptide_ids.inst.get()),
                    deref(convString(experiment_name).get())
                )
            else:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(peptide_ids.inst.get()),
                    deref(convString(experiment_name).get()),
                    deref(lookup.inst.get())
                )
    
    def store(self, filename, protein_ids, peptide_ids, mz_file=None, mz_name=None, 
              peptideprophet_analyzed=False, rt_tolerance=0.0):
        """
        store(self, filename, protein_ids, peptide_ids, mz_file=None, mz_name=None, 
              peptideprophet_analyzed=False, rt_tolerance=0.0)
        
        Stores the identifications to a PepXML file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the PepXML file to store
        :type filename: str or bytes or String
        :param protein_ids: List of protein identifications to store
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container of peptide identifications to store.
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param mz_file: Optional mz file name
        :type mz_file: str or bytes or String, optional
        :param mz_name: Optional mz name
        :type mz_name: str or bytes or String, optional
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
        # Import at function level to avoid circular dependencies
        from . import PeptideIdentificationList as _PeptideIdentificationList
        
        # Check if peptide_ids is a Python list (old API) or PeptideIdentificationList (new API)
        if isinstance(peptide_ids, list):
            # Old API: Convert Python list to PeptideIdentificationList
            temp_peptide_ids = _PeptideIdentificationList()
            temp_peptide_ids.extend(peptide_ids)
            
            if mz_file is None:
                self.inst.get().store(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(temp_peptide_ids.inst.get())
                )
            else:
                self.inst.get().store(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(temp_peptide_ids.inst.get()),
                    deref(convString(mz_file).get()),
                    deref(convString(mz_name).get()),
                    peptideprophet_analyzed,
                    rt_tolerance
                )
        else:
            # New API: Use PeptideIdentificationList directly
            if mz_file is None:
                self.inst.get().store(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(peptide_ids.inst.get())
                )
            else:
                self.inst.get().store(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(peptide_ids.inst.get()),
                    deref(convString(mz_file).get()),
                    deref(convString(mz_name).get()),
                    peptideprophet_analyzed,
                    rt_tolerance
                )
