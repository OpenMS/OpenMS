


    def __str__(self):
        """
        Return a string representation of the Feature object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the Feature object.

        Returns key properties in a readable format similar to:
        Feature(rt=1234.5, mz=445.678, intensity=100000.0, charge=2, quality=0.95)
        """
        cdef double rt = self.getRT()
        cdef double mz = self.getMZ()
        cdef float intensity = self.getIntensity()
        cdef int charge = self.getCharge()
        cdef float quality = self.getOverallQuality()

        # Build the representation string
        parts = []
        parts.append(f"rt={rt:.2f}")
        parts.append(f"mz={mz:.4f}")
        parts.append(f"intensity={intensity:.1f}")
        if charge != 0:
            parts.append(f"charge={charge}")
        # Quality is a core property, always show it (typically in range 0-1)
        parts.append(f"quality={quality:.3f}")

        # Add number of subordinates if present
        subordinates = self.getSubordinates()
        num_subordinates = len(subordinates)
        if num_subordinates > 0:
            parts.append(f"subordinates={num_subordinates}")

        # Add number of peptide IDs if present
        peptide_ids = self.getPeptideIdentifications()
        num_peptide_ids = len(peptide_ids)
        if num_peptide_ids > 0:
            parts.append(f"peptide_ids={num_peptide_ids}")

        return f"Feature({', '.join(parts)})"
