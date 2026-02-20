


    def __str__(self):
        """
        Return a string representation of the ReactionMonitoringTransition object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the ReactionMonitoringTransition object.

        Returns key properties in a readable format:
        ReactionMonitoringTransition(precursor_mz=500.25, product_mz=300.15, peptide_ref='PEPTIDER')
        """
        cdef double precursor_mz = self.getPrecursorMZ()
        cdef double product_mz = self.getProductMZ()
        cdef double lib_intensity = self.getLibraryIntensity()

        # These methods already return Python strings
        native_id_str = self.getNativeID()
        peptide_ref_str = self.getPeptideRef()

        parts = []
        if native_id_str:
            parts.append(f"id='{native_id_str}'")
        parts.append(f"precursor_mz={precursor_mz:.4f}")
        parts.append(f"product_mz={product_mz:.4f}")
        if peptide_ref_str:
            parts.append(f"peptide_ref='{peptide_ref_str}'")
        if lib_intensity > 0:
            parts.append(f"library_intensity={lib_intensity:.1f}")

        return f"ReactionMonitoringTransition({', '.join(parts)})"
