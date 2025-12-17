from UniqueIdInterface cimport setUniqueId as _setUniqueId


    def setUniqueIds(self):
        """
        setUniqueIds(self: ConsensusMap) -> None
        
        Assign unique IDs to all features in the consensus map.
        """
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def __repr__(self):
        """
        __repr__(self: ConsensusMap) -> str
        
        Return a string representation of the ConsensusMap object.

        Returns key properties in a readable format:
        ConsensusMap(num_consensus_features=100, num_maps=3, experiment_type='label-free')
        """
        cdef int num_features = self.inst.get().size()
        
        parts = []
        parts.append(f"num_consensus_features={num_features}")
        
        # Get number of column headers (input maps)
        column_headers = self.getColumnHeaders()
        if len(column_headers) > 0:
            parts.append(f"num_maps={len(column_headers)}")
        
        # Get experiment type
        exp_type = self.getExperimentType()
        if exp_type:
            parts.append(f"experiment_type='{exp_type}'")
        
        # Get number of protein identifications
        protein_ids = self.getProteinIdentifications()
        if len(protein_ids) > 0:
            parts.append(f"protein_ids={len(protein_ids)}")
        
        return f"ConsensusMap({', '.join(parts)})"

    def __len__(self):
        """
        __len__(self: ConsensusMap) -> int

        Return the number of consensus features in the map.

        Enables use of Python's built-in len() function.

        Returns:
            int: The number of consensus features in this map.
        """
        return self.inst.get().size()

    def append(self, ConsensusFeature item):
        """
        append(self: ConsensusMap, item: ConsensusFeature) -> None

        Add a single consensus feature to the end of the map.

        This method provides a Pythonic interface equivalent to push_back().

        Args:
            item: A single ConsensusFeature object to append.
        """
        self.inst.get().push_back(deref(item.inst.get()))

    def extend(self, items):
        """
        extend(self: ConsensusMap, items: Iterable[ConsensusFeature]) -> None

        Add multiple consensus features to the end of the map.

        Args:
            items: Can be:
                - A list/iterable of ConsensusFeature objects
                - Another ConsensusMap object

        Raises:
            TypeError: If items is not iterable or another ConsensusMap.
        """
        if hasattr(items, '__iter__') and not hasattr(items, 'inst'):
            # Handle regular iterables (list, tuple, etc.)
            for consensus_feature in items:
                self.inst.get().push_back(deref((<ConsensusFeature>consensus_feature).inst.get()))
        elif hasattr(items, 'inst') and hasattr(items, '__len__'):
            # Handle another ConsensusMap
            for consensus_feature in items:
                self.inst.get().push_back(deref((<ConsensusFeature>consensus_feature).inst.get()))
        else:
            raise TypeError("extend() argument must be iterable or another ConsensusMap")

    def getColumnHeaders(self):
        """
        getColumnHeaders(self: ConsensusMap) -> Dict[int, ColumnHeader]
        
        Get the column headers describing the input maps.

        # ColumnHeaders is a type alias for Map<..> which can not be
        # handled by autowrap automatically. So we have to provide a manual
        # converter here:
        #
        # (the wrapper works on linux, but msvc complains
        # a lot about the generated code... uwe schmitt)
        """
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
        """
        setColumnHeaders(self: ConsensusMap, in_0: Dict[int, ColumnHeader]) -> None
        
        Set the column headers describing the input maps.
        """
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
