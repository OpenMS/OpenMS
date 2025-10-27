

    def __str__(self):
        """
        Return a string representation of the PeptideHit object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the PeptideHit object.
        
        Returns key properties in a readable format similar to:
        PeptideHit(score=18.1, sequence='PEPTIDER', charge=2, rank=1, protein_refs=['PH_6057'])
        
        If multiple protein evidences exist, also shows aa_before, aa_after, start, end arrays.
        """
        cdef AASequence seq = self.getSequence()
        cdef double score = self.getScore()
        cdef int charge = self.getCharge()
        cdef unsigned int rank = self.getRank()
        
        # Get the sequence as a string
        seq_str = seq.toString().decode('utf-8') if seq.toString() else ""
        
        # Build the representation string
        parts = []
        parts.append(f"score={score}")
        if seq_str:
            parts.append(f"sequence='{seq_str}'")
        parts.append(f"charge={charge}")
        parts.append(f"rank={rank}")
        
        # Add protein evidences if present
        cdef libcpp_vector[PeptideEvidence] evidences = self.getPeptideEvidences()
        if evidences.size() > 0:
            protein_accs = []
            aa_before_list = []
            aa_after_list = []
            start_list = []
            end_list = []
            
            for i in range(evidences.size()):
                acc = evidences[i].getProteinAccession().decode('utf-8')
                if acc:
                    protein_accs.append(acc)
                aa_before_list.append(chr(evidences[i].getAABefore()))
                aa_after_list.append(chr(evidences[i].getAAAfter()))
                start_list.append(evidences[i].getStart())
                end_list.append(evidences[i].getEnd())
            
            # Only include detailed evidence info if we have multiple evidences
            # or if the lists contain useful information
            if protein_accs:
                parts.append(f"protein_refs={protein_accs}")
            
            # Add aa_before, aa_after, start, end if there are evidences
            if evidences.size() > 0:
                # Check if these values are meaningful (not all the same or default)
                if len(set(aa_before_list)) > 1 or (aa_before_list and aa_before_list[0] not in ['X', '\x00']):
                    parts.append(f"aa_before={aa_before_list}")
                if len(set(aa_after_list)) > 1 or (aa_after_list and aa_after_list[0] not in ['X', '\x00']):
                    parts.append(f"aa_after={aa_after_list}")
                if len(set(start_list)) > 1 or (start_list and start_list[0] >= 0):
                    parts.append(f"start={start_list}")
                if len(set(end_list)) > 1 or (end_list and end_list[0] >= 0):
                    parts.append(f"end={end_list}")
        
        return f"PeptideHit({', '.join(parts)})"
