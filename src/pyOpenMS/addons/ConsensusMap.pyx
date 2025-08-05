from UniqueIdInterface cimport setUniqueId as _setUniqueId


    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def getColumnHeaders(self):

        # ColumnHeaders is a type alias for Map<..> which can not be
        # handled by autowrap automatically. So we have to provide a manual
        # converter here:
        #
        # (the wrapper works on linux, but msvc complains
        # a lot about the generated code... uwe schmitt)

        cdef ColumnHeaders _r = self.inst.get().getColumnHeaders()
        py_result = dict()
        cdef ColumnHeaders_iterator it__r = _r.begin()
        cdef ColumnHeader item_py_result
        while it__r != _r.end():
           item_py_result = ColumnHeader.__new__(ColumnHeader)
           item_py_result.inst = shared_ptr[_ColumnHeader](new _ColumnHeader((deref(it__r)).second))
           py_result[<UInt64>(deref(it__r).first)] = item_py_result
           inc(it__r)
        return py_result

    def setColumnHeaders(self, dict in_0 ):
        assert isinstance(in_0, dict) and all(isinstance(k, int) for k in in_0.keys()) and all(isinstance(v, ColumnHeader) for v in in_0.values()), 'arg in_0 wrong type'
        cdef ColumnHeaders v0
        for key, value in in_0.items():
           v0[<UInt64> key] = deref((<ColumnHeader>value).inst.get())

        # we have to utilize AutowrapRefHolder here, because cython does not
        # like
        #
        #     self.inst.get().getColumnHeaders() = v0
        #
        # and if you choose
        #
        #     cdef ColumnHeaders & ref = self.inst.get().getColumnHeaders()
        #     ref = v0
        #
        # the c++ compiler will complain, as cython generates invalid code in
        # this case: Cython creates two statements from the first line:
        #
        #     ColumnHeader & ref;   // invalid: ref not initailized
        #     ref = self.inst.get().getColumnHeaders()
        #
        # wrapping this ref into a class holding the ref and using pointers
        # works. the following code results in c++ code simialar to
        #
        #    AutowrapRefHolder<ColumnHeader> * refholder;
        #    refholder = new AutowrapRefHolder<ColumnHeader>(self.inst...)
        #    refholder.assign(v0)

        cdef AutowrapRefHolder[ColumnHeaders] * refholder
        refholder = new AutowrapRefHolder[ColumnHeaders](self.inst.get().getColumnHeaders())
        refholder.assign(v0)

        cdef replace_in_0 = dict()
        cdef ColumnHeaders_iterator it_in_0 = v0.begin()
        cdef ColumnHeader item_in_0
        while it_in_0 != v0.end():
           item_in_0 = ColumnHeader.__new__(ColumnHeader)
           item_in_0.inst = shared_ptr[_ColumnHeader](new _ColumnHeader((deref(it_in_0)).second))
           replace_in_0[<UInt64> deref(it_in_0).first] = item_in_0
           inc(it_in_0)
        in_0.clear()
        in_0.update(replace_in_0)

    def __len__(self):
        """Return the number of consensus features in the map."""
        return self.inst.get().size()

    def append(self, item):
        """
        Add one or more consensus features to the end of the map.
        
        Args:
            item: Can be:
                - A single ConsensusFeature object
                - A list/iterable of ConsensusFeature objects
                - Another ConsensusMap object
        """
        if hasattr(item, '__iter__') and not hasattr(item, 'inst'):
            # It's a list/iterable of items, but not a pyOpenMS object
            for consensus_feature in item:
                self.inst.get().push_back(deref(consensus_feature.inst.get()))
        elif hasattr(item, 'inst') and hasattr(item, '__len__'):
            # It's another ConsensusMap or similar container object
            for consensus_feature in item:
                self.inst.get().push_back(deref(consensus_feature.inst.get()))
        else:
            # It's a single ConsensusFeature object
            self.inst.get().push_back(deref(item.inst.get()))
