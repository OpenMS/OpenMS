

    def load(self, filename, protein_ids, peptide_ids, *args):
        """
        load(self, filename, protein_ids, peptide_ids, document_id=None)
        
        Loads the identifications from an idXML file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the idXML file to load
        :type filename: str or bytes
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
            # load(filename, protein_ids, peptide_ids, document_id)
            document_id = args[0]
            self.inst.get().load(
                _String(<char*>_bytes(filename)),
                protein_ids,
                cpp_peptide_ids,
                deref(document_id.inst.get())
            )
        else:
            raise ValueError("Invalid number of arguments")
        
        # If we converted from Python list, copy results back
        if needs_conversion:
            peptide_ids[:] = []  # Clear the original list
            for i in range(cpp_peptide_ids.size()):
                peptide_ids.append(cpp_peptide_ids[i])
    
    def store(self, filename, protein_ids, peptide_ids, document_id=b""):
        """
        store(self, filename, protein_ids, peptide_ids, document_id="")
        
        Stores the identifications to an idXML file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param filename: Path to the idXML file to store
        :type filename: str or bytes
        :param protein_ids: List of protein identifications to store
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container of peptide identifications to store.
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param document_id: Optional document identifier
        :type document_id: str or bytes, optional
        
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
        
        # Call the C++ store method
        self.inst.get().store(
            _String(<char*>_bytes(filename)),
            protein_ids,
            cpp_peptide_ids,
            _String(<char*>_bytes(document_id))
        )
