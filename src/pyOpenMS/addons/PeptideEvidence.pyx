

    def __str__(self):
        """
        Return a string representation of the PeptideEvidence object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the PeptideEvidence object.
        
        Returns key properties in a readable format similar to:
        PeptideEvidence(protein='P12345', start=71, end=80, aa_before='R', aa_after='N')
        """
        cdef String acc = self.getProteinAccession()
        cdef int start = self.getStart()
        cdef int end = self.getEnd()
        cdef char aa_before = self.getAABefore()
        cdef char aa_after = self.getAAAfter()
        
        # Get the accession as string
        acc_str = acc.decode('utf-8') if acc.size() > 0 else ""
        
        # Build the representation string
        parts = []
        if acc_str:
            parts.append(f"protein='{acc_str}'")
        if start >= 0:
            parts.append(f"start={start}")
        if end >= 0:
            parts.append(f"end={end}")
        # Only show aa_before/aa_after if meaningful (not null or 'X' for unknown)
        if aa_before != 0 and aa_before != ord('X'):
            parts.append(f"aa_before='{chr(aa_before)}'")
        if aa_after != 0 and aa_after != ord('X'):
            parts.append(f"aa_after='{chr(aa_after)}'")
        
        return f"PeptideEvidence({', '.join(parts)})"
