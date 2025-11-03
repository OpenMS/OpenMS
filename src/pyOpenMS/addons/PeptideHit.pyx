

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
        PeptideHit(score=18.1, sequence='PEPTIDER', charge=2, rank=1, evidences=[...])
        
        Uses PeptideEvidence repr for displaying protein evidence details.
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
        
        # Add protein evidences if present, using their repr
        cdef libcpp_vector[PeptideEvidence] evidences = self.getPeptideEvidences()
        if evidences.size() > 0:
            evidence_reprs = []
            for i in range(evidences.size()):
                # Create Python PeptideEvidence object and use its repr
                py_evidence = PeptideEvidence.__new__(PeptideEvidence)
                py_evidence.inst = shared_ptr[_PeptideEvidence](new _PeptideEvidence(evidences[i]))
                evidence_reprs.append(repr(py_evidence))
            parts.append(f"evidences=[{', '.join(evidence_reprs)}]")
        
        return f"PeptideHit({', '.join(parts)})"
