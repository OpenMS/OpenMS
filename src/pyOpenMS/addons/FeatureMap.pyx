from UniqueIdInterface cimport setUniqueId as _setUniqueId


    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def __len__(self):
        """Return the number of features in the map."""
        return self.inst.get().size()

    def append(self, item):
        """
        Add a single feature to the end of the map.
        
        Args:
            item: A single Feature object
        """
        self.inst.get().push_back(deref(item.inst.get()))

    def extend(self, items):
        """
        Add multiple features to the end of the map.
        
        Args:
            items: Can be:
                - A list/iterable of Feature objects
                - Another FeatureMap object
        """
        if hasattr(items, '__iter__') and not hasattr(items, 'inst'):
            # It's a list/iterable of items, but not a pyOpenMS object
            for feature in items:
                self.inst.get().push_back(deref(feature.inst.get()))
        elif hasattr(items, 'inst') and hasattr(items, '__len__'):
            # It's another FeatureMap or similar container object
            for feature in items:
                self.inst.get().push_back(deref(feature.inst.get()))
        else:
            raise TypeError("extend() argument must be iterable or another FeatureMap")
