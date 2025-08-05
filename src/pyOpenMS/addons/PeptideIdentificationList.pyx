    def __len__(self):
        """Return the number of peptide identifications in the list."""
        return self.inst.get().size()

    def append(self, item):
        """
        Add a single peptide identification to the end of the list.
        
        Args:
            item: A single PeptideIdentification object
        """
        self.inst.get().push_back(deref(item.inst.get()))

    def extend(self, items):
        """
        Add multiple peptide identifications to the end of the list.
        
        Args:
            items: Can be:
                - A list/iterable of PeptideIdentification objects
                - Another PeptideIdentificationList object
        """
        if hasattr(items, '__iter__') and not hasattr(items, 'inst'):
            # It's a list/iterable of items, but not a pyOpenMS object
            for peptide_identification in items:
                self.inst.get().push_back(deref(peptide_identification.inst.get()))
        elif hasattr(items, 'inst') and hasattr(items, '__len__'):
            # It's another PeptideIdentificationList or similar container object
            for peptide_identification in items:
                self.inst.get().push_back(deref(peptide_identification.inst.get()))
        else:
            raise TypeError("extend() argument must be iterable or another PeptideIdentificationList")