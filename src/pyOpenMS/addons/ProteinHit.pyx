

    def __str__(self):
        """
        Return a string representation of the ProteinHit object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the ProteinHit object.

        Returns key properties in a readable format similar to:
        ProteinHit(accession='P12345', score=150.5, coverage=45.2)
        """
        cdef String acc = self.getAccession()
        cdef String desc = self.getDescription()
        cdef double score = self.getScore()
        cdef double coverage = self.getCoverage()

        # Get the accession and description as strings
        acc_str = acc.decode('utf-8') if acc.size() > 0 else ""
        desc_str = desc.decode('utf-8') if desc.size() > 0 else ""

        # Build the representation string
        parts = []
        if acc_str:
            parts.append(f"accession='{acc_str}'")
        parts.append(f"score={score}")
        if coverage >= 0:  # Only include if coverage is meaningful
            parts.append(f"coverage={coverage}")
        if desc_str:
            # Truncate description if too long
            if len(desc_str) > 50:
                desc_str = desc_str[:47] + "..."
            parts.append(f"description='{desc_str}'")
        
        return f"ProteinHit({', '.join(parts)})"
