    def __len__(self):
        """Return the number of peptide identifications in the list."""
        return self.inst.get().size()

    def append(self, item):
        """
        Add one or more peptide identifications to the end of the list.
        
        Args:
            item: Can be:
                - A single PeptideIdentification object
                - A list/iterable of PeptideIdentification objects
                - Another PeptideIdentificationList object
        """
        if hasattr(item, '__iter__') and not hasattr(item, 'inst'):
            # It's a list/iterable of items, but not a pyOpenMS object
            for peptide_identification in item:
                self.inst.get().push_back(deref(peptide_identification.inst.get()))
        elif hasattr(item, 'inst') and hasattr(item, '__len__'):
            # It's another PeptideIdentificationList or similar container object
            for peptide_identification in item:
                self.inst.get().push_back(deref(peptide_identification.inst.get()))
        else:
            # It's a single PeptideIdentification object
            self.inst.get().push_back(deref(item.inst.get()))