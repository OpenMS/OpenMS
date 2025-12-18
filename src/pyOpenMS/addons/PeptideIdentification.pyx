


    def __str__(self):
        """
        Return a string representation of the PeptideIdentification object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the PeptideIdentification object.

        Returns key properties in a readable format similar to:
        PeptideIdentification(rt=1234.5, mz=445.678, score_type='XTandem', num_hits=5)
        """
        # Get score type as string - already returns Python string
        score_type_str = self.getScoreType()

        # Get hits count
        hits = self.getHits()
        num_hits = len(hits)

        # Build the representation string
        parts = []

        # Add RT and MZ if available
        if self.hasRT():
            parts.append(f"rt={self.getRT():.2f}")
        if self.hasMZ():
            parts.append(f"mz={self.getMZ():.4f}")

        if score_type_str:
            parts.append(f"score_type='{score_type_str}'")

        parts.append(f"num_hits={num_hits}")

        # If there are hits, show info about the top hit
        if num_hits > 0:
            top_hit = hits[0]
            top_seq_str = top_hit.getSequence().toString()  # Already returns Python string
            if top_seq_str:
                parts.append(f"top_hit='{top_seq_str}'")
            parts.append(f"top_score={top_hit.getScore()}")

        return f"PeptideIdentification({', '.join(parts)})"
