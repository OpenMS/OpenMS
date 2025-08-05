    def __len__(self):
        """Return the number of peptide identifications in the list."""
        return self.inst.get().size()

    def append(self, peptide_identification):
        """Add a peptide identification to the end of the list."""
        self.inst.get().push_back(peptide_identification.inst.get()[0])