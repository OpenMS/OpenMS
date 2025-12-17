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

    def __len__(self):
        """
        __len__(self: FeatureMap) -> int

        Return the number of features in the map.

        Enables use of Python's built-in len() function.

        Returns:
            int: The number of features in this map.
        """
        return self.inst.get().size()

    def append(self, Feature item):
        """
        append(self: FeatureMap, item: Feature) -> None

        Add a single feature to the end of the map.

        This method provides a Pythonic interface equivalent to push_back().

        Args:
            item: A single Feature object to append.
        """
        self.inst.get().push_back(deref(item.inst.get()))

    def extend(self, items):
        """
        extend(self: FeatureMap, items: Iterable[Feature]) -> None

        Add multiple features to the end of the map.

        Args:
            items: Can be:
                - A list/iterable of Feature objects
                - Another FeatureMap object

        Raises:
            TypeError: If items is not iterable or another FeatureMap.
        """
        if hasattr(items, '__iter__') and not hasattr(items, 'inst'):
            # Handle regular iterables (list, tuple, etc.)
            for feature in items:
                self.inst.get().push_back(deref((<Feature>feature).inst.get()))
        elif hasattr(items, 'inst') and hasattr(items, '__len__'):
            # Handle another FeatureMap
            for feature in items:
                self.inst.get().push_back(deref((<Feature>feature).inst.get()))
        else:
            raise TypeError("extend() argument must be iterable or another FeatureMap")
