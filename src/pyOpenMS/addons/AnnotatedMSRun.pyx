

    def getMapping(AnnotatedMSRun self, size_t idx):
        """Cython signature: Mapping getMapping(size_t idx)"""
        #
        # Get a mapping of spectrum and peptide identification at the given index
        # Parameters
        # ----------
        # idx : size_t
        #   The index of the mapping to retrieve
        #     
        # Returns
        # -------
        # Mapping
        #   A pair containing the spectrum and peptide identification
        #     
        # Raises
        # ------
        # IndexError
        #     If the index is out of bounds        
        assert isinstance(idx, (int, long)), 'arg idx wrong type'
        assert idx < self.size(), 'Requested mapping %s does not exist, there are only %s mappings' % (idx, self.size())
        
        cdef _AnnotatedMSRun * run = self.inst.get()
        cdef _AnnotatedMSRun_Mapping mapping = deref(run)[<size_t>idx]
        
        cdef MSSpectrum py_spectrum = MSSpectrum.__new__(MSSpectrum)
        cdef _MSSpectrum * spec_ptr = new _MSSpectrum(mapping.first)
        py_spectrum.inst = shared_ptr[_MSSpectrum](spec_ptr)
        
        cdef PeptideIdentification py_peptide = PeptideIdentification.__new__(PeptideIdentification)
        cdef _PeptideIdentification * pep_ptr = new _PeptideIdentification(mapping.second)
        py_peptide.inst = shared_ptr[_PeptideIdentification](pep_ptr)
        
        return (py_spectrum, py_peptide)
        
    def __getitem__(AnnotatedMSRun self, size_t idx):
        """Get a mapping of spectrum and peptide identification at the given index"""
        # 
        # Parameters
        # ----------
        # idx : size_t
        #     The index of the mapping to retrieve
        #     
        # Returns
        # -------
        # tuple
        #     A tuple containing (MSSpectrum, PeptideIdentification)
        #     
        # Raises
        # ------
        # IndexError
        #     If the index is out of bounds
        # 
        return self.getMapping(idx)
        
    def __iter__(AnnotatedMSRun self):
        """Iterate over all mappings of spectra and peptide identifications"""
        # 
        # Returns
        # -------
        # iterator
        #     An iterator yielding tuples of (MSSpectrum, PeptideIdentification)        
        cdef size_t n = self.size()
        for i in range(n):
            yield self.getMapping(i)