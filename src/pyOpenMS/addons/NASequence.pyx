




    ## def __iter__(self):
    ##     _r = self.inst.get().getSequence()
    ##     cdef libcpp_vector[const _Ribonucleotide *].iterator it__r = _r.begin()
    ##     cdef Ribonucleotide item_py_result
    ##     while it__r != _r.end():
    ##        item_py_result = Ribonucleotide.__new__(Ribonucleotide)
    ##        item_py_result.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(deref(it__r))))
    ##        yield py_result
    ##        inc(it__r)

    def __iter__(self):

        for r in self.getSequence():
          yield r

    def __len__(self):
        """
        __len__(self: NASequence) -> int
        
        Return the length of the nucleic acid sequence.
        """
        return self.inst.get().size()

    def __str__(self):
        """
        __str__(self: NASequence) -> str
        
        Return the sequence string for human-readable output.
        """
        return self.toString()

    def __repr__(self):
        """
        __repr__(self: NASequence) -> str
        
        Return a string representation of the NASequence object.

        Returns key properties in a readable format:
        NASequence(sequence='ACGU', length=4, mono_mass=1234.56)
        """
        seq_str = self.toString()
        length = self.size()
        mono_mass = self.getMonoWeight()

        parts = []
        if seq_str:
            # Truncate very long sequences (50 chars max, minus 3 for '...')
            max_display_len = 50
            if len(seq_str) > max_display_len:
                parts.append(f"sequence='{seq_str[:max_display_len - 3]}...'")
            else:
                parts.append(f"sequence='{seq_str}'")
        parts.append(f"length={length}")
        parts.append(f"mono_mass={mono_mass:.4f}")

        return f"NASequence({', '.join(parts)})"
