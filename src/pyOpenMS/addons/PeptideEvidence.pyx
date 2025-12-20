

    def __str__(self):
        """
        __str__(self: PeptideEvidence) -> str

        Return a string representation of the PeptideEvidence object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: PeptideEvidence) -> str

        Return a string representation of the PeptideEvidence object.

        Returns key properties in a readable format similar to:
        PeptideEvidence(protein='P12345', start=71, end=80, aa_before='R', aa_after='N')
        """
        # autowrap methods return Python types, not C++ types
        acc_str = self.getProteinAccession()
        start = self.getStart()
        end = self.getEnd()
        aa_before = self.getAABefore()
        aa_after = self.getAAAfter()

        # Build the representation string
        parts = []
        if acc_str:
            parts.append(f"protein='{acc_str}'")
        if start >= 0:
            parts.append(f"start={start}")
        if end >= 0:
            parts.append(f"end={end}")
        # Only show aa_before/aa_after if meaningful (not empty or 'X' for unknown)
        if aa_before and aa_before != 'X':
            parts.append(f"aa_before='{aa_before}'")
        if aa_after and aa_after != 'X':
            parts.append(f"aa_after='{aa_after}'")

        return f"PeptideEvidence({', '.join(parts)})"
