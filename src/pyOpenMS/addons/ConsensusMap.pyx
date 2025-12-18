from UniqueIdInterface cimport setUniqueId as _setUniqueId
cimport numpy as np
import numpy as np
from collections import defaultdict as _defaultdict




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

    def get_df_columns(self, columns='default'):
        """
        get_df_columns(self: ConsensusMap, columns: str = 'default') -> List[str]

        Returns a list of column names that get_df() would produce.

        Useful for discovering available columns before export.

        Args:
            columns (str): 'default' for standard columns, 'all' for all available columns.

        Returns:
            list: List of column name strings.

        Example:
            >>> cmap.get_df_columns()
            ['sequence', 'charge', 'rt', 'mz', 'quality', 'intensity_file1', ...]
        """
        # Metadata columns
        cols = ['sequence', 'charge', 'rt', 'mz', 'quality']

        # Intensity columns depend on experiment type
        labelfree = self.getExperimentType() == "label-free"
        filemeta = self.getColumnHeaders()

        if labelfree:
            # File-wide columns for label-free
            files = list(set([header.filename for header in filemeta.values()]))
            cols.extend(files)
        else:
            # Label columns for labelled experiments
            labels = list(set([header.label for header in filemeta.values()]))
            if len(labels) == 1:
                labels[0] = "intensity"
            cols.extend(labels)
            cols.append('file')

        return cols

    def get_intensity_df(self):
        """
        get_intensity_df(self: ConsensusMap) -> pd.DataFrame

        Generates a pandas DataFrame with feature intensities from each sample in long format (over files).

        For labelled analyses channel intensities will be in one row, therefore resulting in a semi-long/block format.
        Resulting DataFrame can be joined with result from get_metadata_df by their index 'id'.

        Returns:
        pd.DataFrame: intensity DataFrame
        """
        import pandas as pd
        labelfree = self.getExperimentType() == "label-free"
        filemeta = self.getColumnHeaders()  # type: dict[int, ColumnHeader]

        labels = list(set([header.label for header in filemeta.values()]))
        files = list(set([header.filename for header in filemeta.values()]))
        label_to_idx = {k: v for v, k in enumerate(labels)}
        file_to_idx = {k: v for v, k in enumerate(files)}

        def gen(cmap, fun):
            for f in cmap:
                yield from fun(f)

        if not labelfree:

            def extract_row_blocks_channel_wide_file_long(f):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                filerows = _defaultdict(lambda: [0] * len(labels))
                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row = filerows[header.filename]
                    row[label_to_idx[header.label]] = fh.getIntensity()
                return (f.getUniqueId(), filerows)

            def extract_rows_channel_wide_file_long(f):
                uniqueid, rowdict = extract_row_blocks_channel_wide_file_long(f)
                for file, row in rowdict.items():
                    row.append(file)
                    yield tuple([uniqueid] + row)

            if len(labels) == 1:
                labels[0] = "intensity"

            dtypes = [('id', np.dtype('uint64'))] + list(zip(labels, ['f'] * len(labels)))
            dtypes.append(('file', 'U300'))

            # Count actual rows: one row per file per feature (not one per feature)
            total_rows = sum(len(extract_row_blocks_channel_wide_file_long(f)[1]) for f in self)
            intyarr = np.fromiter(iter=gen(self, extract_rows_channel_wide_file_long), dtype=dtypes, count=total_rows)

            return pd.DataFrame(intyarr).set_index('id')

        else:
            # Specialized for LabelFree which has to have only one channel
            def extract_row_blocks_channel_long_file_wide_LF(f):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                row = [0.] * len(files)

                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row[file_to_idx[header.filename]] = fh.getIntensity()

                yield tuple([f.getUniqueId()] + row)

            dtypes = [('id', np.dtype('uint64'))] + list(zip(files, ['f'] * len(files)))

            intyarr = np.fromiter(iter=gen(self, extract_row_blocks_channel_long_file_wide_LF), dtype=dtypes, count=self.size())

            return pd.DataFrame(intyarr).set_index('id')

    def get_metadata_df(self):
        """
        get_metadata_df(self: ConsensusMap) -> pd.DataFrame

        Generates a pandas DataFrame with feature meta data.

        Columns: sequence, charge, rt, mz, quality (indexed by 'id').

        Resulting DataFrame can be joined with result from get_intensity_df by their index 'id'.

        Returns:
            pd.DataFrame: DataFrame with metadata for each feature (sequence, charge,
                         rt, mz, quality). All column names are lowercase snake_case.
        """
        import pandas as pd

        def gen(cmap, fun):
            for f in cmap:
                yield from fun(f)

        def extract_meta_data(f):
            pep = f.getPeptideIdentifications()

            if pep.size() != 0:
                hits = pep[0].getHits()

                if len(hits) != 0:
                    besthit = hits[0]
                    yield f.getUniqueId(), besthit.getSequence().toString(), f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()

                else:
                    yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()

            else:
                yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()

        cnt = self.size()

        mddtypes = [('id', np.dtype('uint64')), ('sequence', 'U200'), ('charge', 'i4'),
                    ('rt', np.dtype('double')), ('mz', np.dtype('double')), ('quality', 'f')]

        mdarr = np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        return pd.DataFrame(mdarr).set_index('id')

    def get_df(self, columns=None):
        """
        get_df(self: ConsensusMap, columns: Optional[List[str]] = None) -> pd.DataFrame

        Generates a pandas DataFrame with both consensus feature meta data and intensities from each sample.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all columns. Use get_df_columns()
                                   to discover available columns.

        Returns:
            pd.DataFrame: meta data and intensity DataFrame

        Example:
            >>> # Get all columns
            >>> df = cmap.get_df()

            >>> # Discover available columns
            >>> print(cmap.get_df_columns())

            >>> # Get only specific columns
            >>> df = cmap.get_df(columns=['sequence', 'mz', 'intensity'])
        """
        import pandas as pd
        if columns is None:
            # No column selection - get everything
            df = pd.concat([self.get_metadata_df(), self.get_intensity_df()], axis=1)
            return df

        # Efficient column selection: only compute what's needed
        requested = set(columns)
        metadata_cols = {'sequence', 'charge', 'rt', 'mz', 'quality'}

        # Get intensity column names
        labelfree = self.getExperimentType() == "label-free"
        filemeta = self.getColumnHeaders()
        if labelfree:
            intensity_cols = set([header.filename for header in filemeta.values()])
        else:
            labels = list(set([header.label for header in filemeta.values()]))
            if len(labels) == 1:
                labels[0] = "intensity"
            intensity_cols = set(labels)
            intensity_cols.add('file')

        need_metadata = len(requested & metadata_cols) > 0
        need_intensity = len(requested & intensity_cols) > 0

        dfs = []
        if need_metadata:
            dfs.append(self.get_metadata_df())
        if need_intensity:
            dfs.append(self.get_intensity_df())

        if not dfs:
            # No columns match - return empty DataFrame with index
            return pd.DataFrame(index=pd.Index([], name='id'))

        if len(dfs) == 1:
            df = dfs[0]
        else:
            df = pd.concat(dfs, axis=1)

        # Filter to requested columns
        available_cols = [c for c in columns if c in df.columns]
        return df[available_cols]

    def to_arrow(self, columns=None):
        """
        to_arrow(self: ConsensusMap, columns: Optional[List[str]] = None) -> pa.Table

        Returns an Apache Arrow Table with consensus feature meta data and intensities.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all columns. Use get_df_columns()
                                   to discover available columns.

        Returns:
            pyarrow.Table: Arrow Table with consensus feature data.

        Example:
            >>> # Get all columns
            >>> table = cmap.to_arrow()

            >>> # Convert to pandas (zero-copy with pandas 2.0+)
            >>> df = table.to_pandas()

            >>> # Convert to polars
            >>> import polars as pl
            >>> df = pl.from_arrow(table)
        """
        import pyarrow as pa
        df = self.get_df(columns=columns)
        return pa.Table.from_pandas(df)
