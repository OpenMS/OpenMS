from UniqueIdInterface cimport setUniqueId as _setUniqueId


    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def __repr__(self):
        """
        __repr__(self: FeatureMap) -> str
        
        Return a string representation of the FeatureMap object.

        Returns key properties in a readable format:
        FeatureMap(num_features=100, protein_ids=2, unassigned_peptide_ids=5)
        """
        cdef int num_features = self.inst.get().size()
        
        parts = []
        parts.append(f"num_features={num_features}")
        
        # Get number of protein identifications
        protein_ids = self.getProteinIdentifications()
        if len(protein_ids) > 0:
            parts.append(f"protein_ids={len(protein_ids)}")
        
        # Get number of unassigned peptide identifications
        unassigned_peptide_ids = self.getUnassignedPeptideIdentifications()
        if len(unassigned_peptide_ids) > 0:
            parts.append(f"unassigned_peptide_ids={len(unassigned_peptide_ids)}")
        
        return f"FeatureMap({', '.join(parts)})"
