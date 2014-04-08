from UniqueIdInterface cimport setUniqueId as _setUniqueId


    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))
