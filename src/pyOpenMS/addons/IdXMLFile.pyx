

    def load(self, filename, protein_ids, peptide_ids, document_id=None):
        """
        load(self, filename, protein_ids, peptide_ids, document_id=None)
        
        Loads the identifications from an idXML file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the idXML file to load
        :type filename: str or bytes or String
        :param protein_ids: List to store protein identifications (modified in place)
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param document_id: Optional document identifier (output parameter)
        :type document_id: String, optional
        
        Example::
        
            # New style (3.5+)
            protein_ids = []
            peptide_ids = pyopenms.PeptideIdentificationList()
            pyopenms.IdXMLFile().load("test.idXML", protein_ids, peptide_ids)
            
            # Old style (backward compatible)
            protein_ids = []
            peptide_ids = []
            pyopenms.IdXMLFile().load("test.idXML", protein_ids, peptide_ids)
        """
        # Import at function level to avoid circular dependencies
        from . import PeptideIdentificationList as _PeptideIdentificationList
        
        # Check if peptide_ids is a Python list (old API) or PeptideIdentificationList (new API)
        if isinstance(peptide_ids, list):
            # Old API: Convert Python list to PeptideIdentificationList
            temp_peptide_ids = _PeptideIdentificationList()
            temp_peptide_ids.extend(peptide_ids)
            
            # Call the C++ load method
            if document_id is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(temp_peptide_ids.inst.get())
                )
            else:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(temp_peptide_ids.inst.get()),
                    deref(convString(document_id).get())
                )
            
            # Copy results back to Python list
            peptide_ids[:] = list(temp_peptide_ids)
        else:
            # New API: Use PeptideIdentificationList directly
            if document_id is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(peptide_ids.inst.get())
                )
            else:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    protein_ids,
                    deref(peptide_ids.inst.get()),
                    deref(convString(document_id).get())
                )
    
    def store(self, filename, protein_ids, peptide_ids, document_id=b""):
        """
        store(self, filename, protein_ids, peptide_ids, document_id="")
        
        Stores the identifications to an idXML file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the idXML file to store
        :type filename: str or bytes or String
        :param protein_ids: List of protein identifications to store
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container of peptide identifications to store.
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param document_id: Optional document identifier
        :type document_id: str or bytes or String, optional
        
        Example::
        
            # New style (3.5+)
            protein_ids = [...]
            peptide_ids = pyopenms.PeptideIdentificationList()
            peptide_ids.extend([...])
            pyopenms.IdXMLFile().store("test.idXML", protein_ids, peptide_ids)
            
            # Old style (backward compatible)
            protein_ids = [...]
            peptide_ids = [...]
            pyopenms.IdXMLFile().store("test.idXML", protein_ids, peptide_ids)
        """
        # Import at function level to avoid circular dependencies
        from . import PeptideIdentificationList as _PeptideIdentificationList
        
        # Check if peptide_ids is a Python list (old API) or PeptideIdentificationList (new API)
        if isinstance(peptide_ids, list):
            # Old API: Convert Python list to PeptideIdentificationList
            temp_peptide_ids = _PeptideIdentificationList()
            temp_peptide_ids.extend(peptide_ids)
            
            self.inst.get().store(
                deref(convString(filename).get()),
                protein_ids,
                deref(temp_peptide_ids.inst.get()),
                deref(convString(document_id).get())
            )
        else:
            # New API: Use PeptideIdentificationList directly
            self.inst.get().store(
                deref(convString(filename).get()),
                protein_ids,
                deref(peptide_ids.inst.get()),
                deref(convString(document_id).get())
            )
