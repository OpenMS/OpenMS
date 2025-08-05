from UniqueIdInterface cimport setUniqueId as _setUniqueId


    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def __len__(self):
        """Return the number of features in the map."""
        return self.inst.get().size()

    def append(self, feature):
        """Add a feature to the end of the map."""
        self.inst.get().push_back(feature.inst.get()[0])
