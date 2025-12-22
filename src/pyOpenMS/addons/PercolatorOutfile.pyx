

    def load(self, filename, proteins, peptides, lookup, output_score):
        """
        load(self, filename, proteins, peptides, lookup, output_score)
        
        Loads a Percolator output file.
        
        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptides parameter.
        
        :param filename: Path to the Percolator output file to load
        :type filename: str or bytes
        :param proteins: Protein identification object to store results (modified in place)
        :type proteins: ProteinIdentification
        :param peptides: Container to store peptide identifications (modified in place).
                        Accepts both list[PeptideIdentification] (pre-3.5 style) and
                        PeptideIdentificationList (3.5+ style)
        :type peptides: list[PeptideIdentification] or PeptideIdentificationList
        :param lookup: Spectrum metadata lookup for RT and m/z information
        :type lookup: SpectrumMetaDataLookup
        :param output_score: Score type to use for output
        :type output_score: PercolatorOutfile.ScoreType
        
        Example::
        
            # New style (3.5+)
            proteins = pyopenms.ProteinIdentification()
            peptides = pyopenms.PeptideIdentificationList()
            lookup = pyopenms.SpectrumMetaDataLookup()
            score_type = pyopenms.PercolatorOutfile.ScoreType.QVALUE
            pyopenms.PercolatorOutfile().load("test.psms", proteins, peptides, lookup, score_type)
            
            # Old style (backward compatible)
            proteins = pyopenms.ProteinIdentification()
            peptides = []
            lookup = pyopenms.SpectrumMetaDataLookup()
            score_type = pyopenms.PercolatorOutfile.ScoreType.QVALUE
            pyopenms.PercolatorOutfile().load("test.psms", proteins, peptides, lookup, score_type)
        """
        cdef PeptideIdentificationList cpp_peptides
        cdef bool needs_conversion = False
        
        # Check if peptides is a Python list (old API) or PeptideIdentificationList (new API)
        if isinstance(peptides, list):
            # Old API: Convert Python list to PeptideIdentificationList
            needs_conversion = True
            cpp_peptides = PeptideIdentificationList()
            # Copy existing items from Python list to PeptideIdentificationList
            for pep_id in peptides:
                cpp_peptides.push_back(pep_id)
        else:
            # New API: Use PeptideIdentificationList directly
            cpp_peptides = peptides
        
        # Call the C++ load method
        self.inst.get().load(
            _String(<char*>_bytes(filename)),
            deref(proteins.inst.get()),
            cpp_peptides,
            deref(lookup.inst.get()),
            <PercolatorOutfile_ScoreType>output_score
        )
        
        # If we converted from Python list, copy results back
        if needs_conversion:
            peptides[:] = []  # Clear the original list
            for i in range(cpp_peptides.size()):
                peptides.append(cpp_peptides[i])
