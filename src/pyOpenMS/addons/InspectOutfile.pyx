

    def load(self, result_filename, peptide_ids, protein_identification, p_value_threshold, database_filename):
        """
        load(self, result_filename, peptide_ids, protein_identification, p_value_threshold, database_filename)
        
        Loads the results of an Inspect search.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.
        
        :param result_filename: Input file name
        :type result_filename: str or bytes
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param protein_identification: Protein identification object to store results (modified in place)
        :type protein_identification: ProteinIdentification
        :param p_value_threshold: P-value threshold for filtering
        :type p_value_threshold: float
        :param database_filename: Database file name
        :type database_filename: str or bytes
        :return: Vector of record indices
        :rtype: list[int]
        
        Example::
        
            # New style (3.5+)
            peptide_ids = pyopenms.PeptideIdentificationList()
            protein_id = pyopenms.ProteinIdentification()
            records = pyopenms.InspectOutfile().load("test.txt", peptide_ids, protein_id, 0.01, "db.fasta")
            
            # Old style (backward compatible)
            peptide_ids = []
            protein_id = pyopenms.ProteinIdentification()
            records = pyopenms.InspectOutfile().load("test.txt", peptide_ids, protein_id, 0.01, "db.fasta")
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
        cdef libcpp_vector[size_t] result = self.inst.get().load(
            _String(<char*>_bytes(result_filename)),
            cpp_peptide_ids,
            deref(protein_identification.inst.get()),
            p_value_threshold,
            _String(<char*>_bytes(database_filename))
        )
        
        # If we converted from Python list, copy results back
        if needs_conversion:
            peptide_ids[:] = []  # Clear the original list
            for i in range(cpp_peptide_ids.size()):
                peptide_ids.append(cpp_peptide_ids[i])
        
        # Convert result vector to Python list
        return [result[i] for i in range(result.size())]
