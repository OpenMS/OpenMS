from UniqueIdInterface cimport setUniqueId as _setUniqueId


    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def __len__(self):
        """Return the number of features in the map."""
        return self.inst.get().size()

    def append(self, item):
        """
        Add one or more features to the end of the map.
        
        Args:
            item: Can be:
                - A single Feature object
                - A list/iterable of Feature objects
                - Another FeatureMap object
        """
        if hasattr(item, '__iter__') and not hasattr(item, 'inst'):
            # It's a list/iterable of items, but not a pyOpenMS object
            for feature in item:
                self.inst.get().push_back(deref(feature.inst.get()))
        elif hasattr(item, 'inst') and hasattr(item, '__len__'):
            # It's another FeatureMap or similar container object
            for feature in item:
                self.inst.get().push_back(deref(feature.inst.get()))
        else:
            # It's a single Feature object
            self.inst.get().push_back(deref(item.inst.get()))
