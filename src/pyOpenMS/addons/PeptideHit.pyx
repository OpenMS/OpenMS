


    def __str__(self):
        """
        __str__(self: PeptideHit) -> str

        Return a string representation of the PeptideHit object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: PeptideHit) -> str

        Return a string representation of the PeptideHit object.

        Returns key properties in a readable format similar to:
        PeptideHit(score=18.1, sequence='PEPTIDER', charge=2, evidences=[...])
        """
        score = self.getScore()
        charge = self.getCharge()

        # Get the sequence as a string using Python method
        seq = self.getSequence()
        seq_str = seq.toString()  # Already returns a Python string

        # Build the representation string
        parts = []
        parts.append(f"score={score}")
        if seq_str:
            parts.append(f"sequence='{seq_str}'")
        parts.append(f"charge={charge}")

        # Add protein evidences with details
        evidences = self.getPeptideEvidences()
        if len(evidences) > 0:
            ev_strs = [repr(e) for e in evidences]
            parts.append(f"evidences=[{', '.join(ev_strs)}]")

        return f"PeptideHit({', '.join(parts)})"
