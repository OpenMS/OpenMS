from libcpp.map cimport map as libcpp_map


    def getAAFrequencies(self, dict mmap):
        """
        getAAFrequencies(self: AASequence, mmap: dict) -> None
        
        Get amino acid frequencies and populate the provided dictionary.
        """
        cdef libcpp_map[_String, size_t] c_mmap
        self.inst.get().getAAFrequencies(c_mmap)
        for k,v in mmap.iteritems():
            v = c_mmap[ _String(<char *>k) ]

    def __iter__(self):
        """
        __iter__(self: AASequence) -> Iterator[Residue]
        
        Iterate over residues in the amino acid sequence.
        """
        cdef unsigned int n = self.inst.get().size()
        cdef unsigned int i = 0

        cdef Residue py_result
        while i < n:

            py_result = Residue.__new__(Residue)
            py_result.inst = shared_ptr[_Residue](new _Residue( deref(self.inst.get())[i] ))

            yield py_result

            i += 1

    def __len__(self):
        """
        __len__(self: AASequence) -> int
        
        Return the length of the amino acid sequence.
        """
        return self.inst.get().size()

    def __str__(self):
        """
        __str__(self: AASequence) -> str
        
        Return a string representation of the AASequence object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: AASequence) -> str
        
        Return a string representation of the AASequence object.

        Returns key properties in a readable format:
        AASequence(sequence='PEPTM(Oxidation)IDE', length=9, mono_mass=1050.42, modified=True)
        """
        # Use Python methods to get values
        seq_decoded = self.toString()  # Already returns a Python string
        length = self.size()
        mono_mass = self.getMonoWeight()
        is_modified = self.isModified()

        parts = []
        if seq_decoded:
            # Truncate very long sequences
            if len(seq_decoded) > 50:
                parts.append(f"sequence='{seq_decoded[:47]}...'")
            else:
                parts.append(f"sequence='{seq_decoded}'")
        parts.append(f"length={length}")
        parts.append(f"mono_mass={mono_mass:.4f}")
        if is_modified:
            parts.append("modified=True")

        return f"AASequence({', '.join(parts)})"

