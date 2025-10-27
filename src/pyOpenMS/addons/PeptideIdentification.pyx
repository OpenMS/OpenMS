

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
        cdef libcpp_vector[PeptideHit] hits = self.getHits()
        cdef String score_type = self.getScoreType()
        
        # Get score type as string
        score_type_str = score_type.decode('utf-8') if score_type.size() > 0 else ""
        
        # Build the representation string
        parts = []
        
        # Add RT and MZ if available
        if self.hasRT():
            parts.append(f"rt={self.getRT():.2f}")
        if self.hasMZ():
            parts.append(f"mz={self.getMZ():.4f}")
        
        if score_type_str:
            parts.append(f"score_type='{score_type_str}'")
        
        parts.append(f"num_hits={hits.size()}")
        
        # If there are hits, show info about the top hit
        if hits.size() > 0:
            cdef AASequence top_seq = hits[0].getSequence()
            top_seq_str = top_seq.toString().decode('utf-8') if top_seq.toString() else ""
            if top_seq_str:
                parts.append(f"top_hit='{top_seq_str}'")
            parts.append(f"top_score={hits[0].getScore()}")
        
        return f"PeptideIdentification({', '.join(parts)})"
