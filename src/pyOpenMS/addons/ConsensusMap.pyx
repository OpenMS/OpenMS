from UniqueIdInterface cimport setUniqueId as _setUniqueId
cimport numpy as np
import numpy as np
import warnings
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

        :return: The number of consensus features in this map.
        :rtype: int
        """
        return self.inst.get().size()

    def findProteinIdentification(self, identifier):
        """
        Find a ProteinIdentification by its identifier.

        :param identifier: The identifier string to search for
        :return: The matching ProteinIdentification, or None if not found
        :rtype: ProteinIdentification or None
        """
        for prot_id in self.getProteinIdentifications():
            if prot_id.getIdentifier() == identifier:
                return prot_id
        return None

    def append(self, ConsensusFeature item):
        """
        append(self: ConsensusMap, item: ConsensusFeature) -> None

        Add a single consensus feature to the end of the map.

        This method provides a Pythonic interface equivalent to push_back().

        :param item: A single ConsensusFeature object to append.
        :type item: ConsensusFeature
        """
        self.inst.get().push_back(deref(item.inst.get()))

    def extend(self, items):
        """
        extend(self: ConsensusMap, items: Iterable[ConsensusFeature]) -> None

        Add multiple consensus features to the end of the map.

        :param items: Can be:
                      - A list/iterable of ConsensusFeature objects
                      - Another ConsensusMap object
        :type items: Iterable[ConsensusFeature]
        :raises TypeError: If items is not iterable or another ConsensusMap.
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

    def df_columns(self, columns='default'):
        """
        df_columns(self: ConsensusMap, columns: str = 'default') -> List[str]

        Returns a list of column names that to_df() would produce.

        Useful for discovering available columns before export.

        :param columns: 'default' for standard columns, 'all' for all available columns.
        :type columns: str
        :return: List of column name strings.
        :rtype: list

        Example::

            >>> cmap.df_columns()
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

        :return: Intensity DataFrame
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for get_intensity_df(). "
                "Please install it with: pip install pandas"
            )
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

        :return: DataFrame with metadata for each feature (sequence, charge,
                 rt, mz, quality). All column names are lowercase snake_case.
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for get_metadata_df(). "
                "Please install it with: pip install pandas"
            )

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

    def to_df(self, columns=None):
        """
        to_df(self: ConsensusMap, columns: Optional[List[str]] = None) -> pd.DataFrame

        Generates a pandas DataFrame with both consensus feature meta data and intensities from each sample.

        :param columns: List of column names to include. If None,
                        includes all columns. Use df_columns()
                        to discover available columns.
        :type columns: Optional[List[str]]

        :return: Meta data and intensity DataFrame
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed

        Example::

            # Get all columns
            df = consensusmap.to_df()

            # Discover available columns
            print(consensusmap.df_columns())

            # Get only specific columns
            df = consensusmap.to_df(columns=['sequence', 'mz', 'intensity'])
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for to_df(). "
                "Please install it with: pip install pandas"
            )
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

        :param columns: List of column names to include. If None,
                        includes all columns. Use df_columns()
                        to discover available columns.
        :type columns: Optional[List[str]]

        :return: Arrow Table with consensus feature data.
        :rtype: pyarrow.Table

        :raises ImportError: If pyarrow is not installed

        Example::

            # Get all columns
            table = consensusmap.to_arrow()

            # Convert to pandas (zero-copy with pandas 2.0+)
            df = table.to_pandas()

            # Convert to polars
            import polars as pl
            df = pl.from_arrow(table)
        """
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_arrow(). "
                "Please install it with: pip install pyarrow"
            )
        df = self.to_df(columns=columns)
        return pa.Table.from_pandas(df)

    # Deprecated aliases for backward compatibility with pyopenms 3.5.0
    def get_df(self, *args, **kwargs):
        """Deprecated: Use to_df() instead."""
        import warnings
        warnings.warn(
            "get_df() is deprecated and will be removed in a future version. "
            "Use to_df() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.to_df(*args, **kwargs)

    def get_df_columns(self, *args, **kwargs):
        """Deprecated: Use df_columns() instead."""
        import warnings
        warnings.warn(
            "get_df_columns() is deprecated and will be removed in a future version. "
            "Use df_columns() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        return self.df_columns(*args, **kwargs)

    def to_feature_arrow(self, reference_file_name="", columns=None,
                          include_modifications=True, scan_format="scan"):
        """
        to_feature_arrow(self: ConsensusMap, reference_file_name: str = "", columns: list = None, include_modifications: bool = True, scan_format: str = "scan") -> pa.Table

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export consensus features as Apache Arrow Table following the QPX feature schema.

        Each ConsensusFeature becomes one row. The best PeptideHit (rank 0) is used
        for identification columns. Intensities from FeatureHandles are mapped to
        structured intensity lists using ColumnHeaders for sample/channel mapping.

        Results are sorted by rt, observed_mz, precursor_charge for consistent ordering.

        :param reference_file_name: Source file name to include in each row.
        :type reference_file_name: str
        :param columns: List of column names to include. If None, includes all.
                        Use feature_columns() to discover available columns.
        :type columns: list
        :param include_modifications: Include modification details. Default True.
        :type include_modifications: bool
        :param scan_format: Controls the scan column: "scan" extracts scan number,
                            "nativeId" uses raw native ID string.
        :type scan_format: str
        :return: Arrow Table with QPX feature schema.
        :rtype: pyarrow.Table
        :raises ImportError: If pyarrow is not installed
        """
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_feature_arrow(). "
                "Please install it with: pip install pyarrow"
            )

        if scan_format not in ("scan", "nativeId"):
            raise ValueError(f"scan_format must be 'scan' or 'nativeId', got '{scan_format}'")

        from . import SpectrumLookup as _SpectrumLookup
        from . import IDScoreSwitcherAlgorithm as _IDScoreSwitcherAlgorithm
        from . import Scores as _Scores
        from . import ProForma as _ProForma
        _IDType = _Scores.IDType

        # Native ID type accessions for scan number extraction
        _native_id_accessions = [
            "MS:1000768", "MS:1000769", "MS:1000771", "MS:1001508",
            "MS:1000774", "MS:1000777", "MS:1001530",
        ]
        _spectrum_lookup = _SpectrumLookup()

        def _extract_scan_number(native_id):
            if not native_id:
                return None
            for accession in _native_id_accessions:
                try:
                    scan_num = _spectrum_lookup.extractScanNumber(native_id, accession)
                    if scan_num >= 0:
                        return scan_num
                except Exception:
                    continue
            return None

        # Build column header lookup: map_index -> (filename, label)
        filemeta = self.getColumnHeaders()

        # Build protein group q-value lookup from ProteinIdentifications
        pg_qvalue_lookup = {}  # accession -> q-value
        gg_accessions_lookup = {}  # accession -> list of gene accessions
        gg_names_lookup = {}  # accession -> list of gene names

        for prot_id in self.getProteinIdentifications():
            # Build protein hit lookup for gene info
            for ph in prot_id.getHits():
                acc = ph.getAccession()
                # Gene info from metavalues if available
                gene_accessions = []
                gene_names = []
                if ph.metaValueExists("gene_accession"):
                    val = ph.getMetaValue("gene_accession")
                    gene_accessions.append(str(val))
                if ph.metaValueExists("gene_name"):
                    val = ph.getMetaValue("gene_name")
                    gene_names.append(str(val))
                if gene_accessions:
                    gg_accessions_lookup[acc] = gene_accessions
                if gene_names:
                    gg_names_lookup[acc] = gene_names

            # Protein groups with q-values
            for pg in prot_id.getProteinGroups():
                qvalue = pg.probability
                # pg.accessions are bytes, decode to str for lookup
                accessions = [a.decode('utf-8') if isinstance(a, bytes) else a for a in pg.accessions]
                for acc in accessions:
                    pg_qvalue_lookup[acc] = qvalue

            # Indistinguishable proteins
            for ig in prot_id.getIndistinguishableProteins():
                # ig.accessions are bytes, decode to str for lookup
                accessions = [a.decode('utf-8') if isinstance(a, bytes) else a for a in ig.accessions]
                qvalue = ig.probability
                for acc in accessions:
                    if acc not in pg_qvalue_lookup:
                        pg_qvalue_lookup[acc] = qvalue

        # Create IDScoreSwitcherAlgorithm instance
        idsa = _IDScoreSwitcherAlgorithm()

        # Initialize column lists
        all_sequence = []
        all_peptidoform = []
        all_modifications = []
        all_precursor_charge = []
        all_calculated_mz = []
        all_observed_mz = []
        all_rt = []
        all_pep = []
        all_is_decoy = []
        all_additional_scores = []
        all_predicted_rt = []
        all_reference_file = []
        all_cv_params = []
        all_scan = []
        all_ion_mobility = []
        all_start_ion_mobility = []
        all_stop_ion_mobility = []
        all_intensities = []
        all_additional_intensities = []
        all_pg_accessions = []
        all_anchor_protein = []
        all_unique = []
        all_pg_global_qvalue = []
        all_gg_accessions = []
        all_gg_names = []
        all_scan_reference_file_name = []
        all_rt_start = []
        all_rt_stop = []
        # OpenMS-specific
        all_quality = []
        all_score = []
        all_score_type = []
        all_spectrum_reference = []
        all_feature_metavalues = []

        # Excluded metavalues (have dedicated columns)
        _excluded_metavalues = {"target_decoy", "predicted_RT", "predicted_rt"}

        def _get_value_type(val):
            if type(val).__name__ == 'bool':
                return "bool"
            elif isinstance(val, int):
                return "int"
            elif isinstance(val, float):
                return "float"
            else:
                return "str"

        for cf in self:
            pep_ids = cf.getPeptideIdentifications()
            has_id = len(pep_ids) > 0
            best_hit = None
            seq = None
            charge = cf.getCharge()

            if has_id:
                hits = pep_ids[0].getHits()
                if len(hits) > 0:
                    best_hit = hits[0]
                    seq = best_hit.getSequence()

            # sequence / peptidoform (ProForma notation)
            if seq is not None:
                all_sequence.append(seq.toUnmodifiedString())
                pf = _ProForma.fromAASequence(seq)
                all_peptidoform.append(_ProForma.toString(pf, _ProForma.WriteMode.CANONICAL))
            else:
                all_sequence.append(None)
                all_peptidoform.append(None)

            # modifications (duplicated extraction pattern from PSM export)
            modifications = []
            if include_modifications and seq is not None and seq.isModified():
                mod_dict = {}

                if seq.hasNTerminalModification():
                    mod = seq.getNTerminalModification()
                    mod_name = seq.getNTerminalModificationName()
                    accession = mod.getUniModRecordId() if mod else None
                    accession_str = f"UNIMOD:{accession}" if accession and accession > 0 else None
                    position_str = "N-term.0"
                    if mod_name not in mod_dict:
                        mod_dict[mod_name] = {"name": mod_name, "accession": accession_str, "positions": []}
                    mod_dict[mod_name]["positions"].append({"position": position_str, "scores": None})

                for pos, residue in enumerate(seq):
                    if residue.isModified():
                        res_mod = residue.getModification()
                        mod_name = residue.getModificationName()
                        accession = res_mod.getUniModRecordId() if res_mod else None
                        accession_str = f"UNIMOD:{accession}" if accession and accession > 0 else None
                        aa_code = residue.getOneLetterCode()
                        position_str = f"{aa_code}.{pos + 1}"
                        if mod_name not in mod_dict:
                            mod_dict[mod_name] = {"name": mod_name, "accession": accession_str, "positions": []}
                        mod_dict[mod_name]["positions"].append({"position": position_str, "scores": None})

                if seq.hasCTerminalModification():
                    mod = seq.getCTerminalModification()
                    mod_name = seq.getCTerminalModificationName()
                    accession = mod.getUniModRecordId() if mod else None
                    accession_str = f"UNIMOD:{accession}" if accession and accession > 0 else None
                    position_str = f"C-term.{seq.size() + 1}"
                    if mod_name not in mod_dict:
                        mod_dict[mod_name] = {"name": mod_name, "accession": accession_str, "positions": []}
                    mod_dict[mod_name]["positions"].append({"position": position_str, "scores": None})

                modifications = list(mod_dict.values())
            all_modifications.append(modifications)

            # precursor_charge
            all_precursor_charge.append(charge)

            # calculated_mz
            calculated_mz = None
            if seq is not None and charge > 0:
                try:
                    calculated_mz = seq.getMZ(charge)
                except Exception:
                    pass
            all_calculated_mz.append(calculated_mz)

            # observed_mz, rt
            all_observed_mz.append(cf.getMZ())
            all_rt.append(cf.getRT())

            # PEP
            pep_value = None
            if has_id:
                pep_search_result = idsa.findScoreType(pep_ids[0], _IDType.PEP)
                pep_score_name = pep_search_result.score_name
                if pep_score_name and best_hit is not None:
                    if pep_search_result.is_main_score_type:
                        pep_value = best_hit.getScore()
                    elif best_hit.metaValueExists(pep_score_name):
                        pep_value = best_hit.getMetaValue(pep_score_name)
            all_pep.append(pep_value)

            # is_decoy
            is_decoy = None
            if best_hit is not None and best_hit.metaValueExists("target_decoy"):
                td = best_hit.getMetaValue("target_decoy")
                is_decoy = 1 if str(td).startswith("decoy") else 0
            all_is_decoy.append(is_decoy)

            # additional_scores from best hit
            additional_scores = []
            if best_hit is not None:
                keys = []
                best_hit.getKeys(keys)
                for key in keys:
                    # Normalize key to string (getKeys may return bytes in some cases)
                    key_str = key.decode("utf-8") if isinstance(key, bytes) else key
                    if key_str not in _excluded_metavalues:
                        val = best_hit.getMetaValue(key)
                        if isinstance(val, (int, float)) and _Scores.isKnownScoreType(key_str):
                            higher_better = None
                            try:
                                score_type_enum = _IDScoreSwitcherAlgorithm.toScoreTypeEnum(key_str)
                                higher_better = idsa.isScoreTypeHigherBetter(score_type_enum)
                            except Exception:
                                pass
                            additional_scores.append({
                                "score_name": key_str,
                                "score_value": float(val),
                                "higher_better": higher_better
                            })
            all_additional_scores.append(additional_scores)

            # predicted_rt
            predicted_rt = None
            if best_hit is not None:
                if best_hit.metaValueExists("predicted_RT"):
                    predicted_rt = best_hit.getMetaValue("predicted_RT")
                elif best_hit.metaValueExists("predicted_rt"):
                    predicted_rt = best_hit.getMetaValue("predicted_rt")
            all_predicted_rt.append(predicted_rt)

            # reference_file_name
            all_reference_file.append(reference_file_name)

            # cv_params (reserved)
            all_cv_params.append(None)

            # scan + scan_reference_file_name from best PSM's spectrum reference
            spec_ref = ""
            if has_id and pep_ids[0].metaValueExists("spectrum_reference"):
                spec_ref = pep_ids[0].getMetaValue("spectrum_reference")

            if scan_format == "nativeId":
                all_scan.append(spec_ref if spec_ref else None)
            else:
                scan_num = _extract_scan_number(spec_ref)
                all_scan.append(str(scan_num) if scan_num is not None else None)

            all_scan_reference_file_name.append(reference_file_name if spec_ref else None)

            # ion_mobility
            ion_mobility = None
            if has_id:
                if pep_ids[0].metaValueExists("ion_mobility"):
                    ion_mobility = pep_ids[0].getMetaValue("ion_mobility")
                elif pep_ids[0].metaValueExists("IM"):
                    ion_mobility = pep_ids[0].getMetaValue("IM")
            all_ion_mobility.append(ion_mobility)

            # start/stop ion mobility from metavalues
            start_im = None
            stop_im = None
            if cf.metaValueExists("start_ion_mobility"):
                start_im = cf.getMetaValue("start_ion_mobility")
            if cf.metaValueExists("stop_ion_mobility"):
                stop_im = cf.getMetaValue("stop_ion_mobility")
            all_start_ion_mobility.append(start_im)
            all_stop_ion_mobility.append(stop_im)

            # intensities from FeatureHandles + ColumnHeaders
            intensities = []
            for fh in cf.getFeatureList():
                map_idx = fh.getMapIndex()
                if map_idx in filemeta:
                    header = filemeta[map_idx]
                    intensities.append({
                        "sample_accession": header.filename,
                        "channel": header.label,
                        "intensity": fh.getIntensity()
                    })
            all_intensities.append(intensities)
            all_additional_intensities.append(None)

            # protein info from best hit - deduplicate while preserving order
            protein_accessions = []
            seen_accs = set()
            if best_hit is not None:
                evidences = best_hit.getPeptideEvidences()
                for ev in evidences:
                    acc = ev.getProteinAccession()
                    if acc not in seen_accs:
                        seen_accs.add(acc)
                        protein_accessions.append(acc)

            all_pg_accessions.append(protein_accessions)
            all_anchor_protein.append(protein_accessions[0] if protein_accessions else None)
            all_unique.append(1 if len(protein_accessions) == 1 else 0)

            # pg_global_qvalue - use minimum q-value across all protein accessions
            pg_qval = None
            if protein_accessions:
                min_qval = None
                for acc in protein_accessions:
                    if acc in pg_qvalue_lookup:
                        qval = pg_qvalue_lookup[acc]
                        if min_qval is None or qval < min_qval:
                            min_qval = qval
                pg_qval = min_qval
            all_pg_global_qvalue.append(pg_qval)

            # gene info
            gg_acc = []
            gg_nam = []
            for acc in protein_accessions:
                if acc in gg_accessions_lookup:
                    gg_acc.extend(gg_accessions_lookup[acc])
                if acc in gg_names_lookup:
                    gg_nam.extend(gg_names_lookup[acc])
            all_gg_accessions.append(list(set(gg_acc)) if gg_acc else [])
            all_gg_names.append(list(set(gg_nam)) if gg_nam else [])

            # rt_start / rt_stop from feature width or metavalues
            rt_start = None
            rt_stop = None
            if cf.metaValueExists("rt_start"):
                rt_start = cf.getMetaValue("rt_start")
            if cf.metaValueExists("rt_stop"):
                rt_stop = cf.getMetaValue("rt_stop")
            if rt_start is None and rt_stop is None:
                width = cf.getWidth()
                if width > 0:
                    rt_start = cf.getRT() - width / 2.0
                    rt_stop = cf.getRT() + width / 2.0
            all_rt_start.append(rt_start)
            all_rt_stop.append(rt_stop)

            # OpenMS-specific columns
            all_quality.append(cf.getQuality())

            score_val = None
            score_type_str = None
            if has_id and best_hit is not None:
                score_val = best_hit.getScore()
                score_type_str = pep_ids[0].getScoreType()
            all_score.append(score_val)
            all_score_type.append(score_type_str)

            all_spectrum_reference.append(spec_ref if spec_ref else None)

            # feature_metavalues
            feature_mvs = []
            cf_keys = []
            cf.getKeys(cf_keys)
            _excluded_feature_mvs = {"rt_start", "rt_stop", "start_ion_mobility", "stop_ion_mobility"}
            for key in cf_keys:
                # Normalize key to string (getKeys may return bytes in some cases)
                key_str = key.decode("utf-8") if isinstance(key, bytes) else key
                if key_str not in _excluded_feature_mvs:
                    val = cf.getMetaValue(key)
                    val_type = _get_value_type(val)
                    feature_mvs.append({"name": key_str, "value": val, "value_type": val_type})
            all_feature_metavalues.append(feature_mvs)

        # Build Arrow table
        columns_set = set(columns) if columns is not None else None

        def should_include(col_name):
            return columns_set is None or col_name in columns_set

        # Define nested types
        position_type = pa.struct([
            ("position", pa.utf8()),
            ("scores", pa.float64())
        ])
        modification_type = pa.struct([
            ("name", pa.utf8()),
            ("accession", pa.utf8()),
            ("positions", pa.list_(position_type))
        ])
        score_type_struct = pa.struct([
            ("score_name", pa.utf8()),
            ("score_value", pa.float64()),
            ("higher_better", pa.bool_())
        ])
        intensity_type = pa.struct([
            ("sample_accession", pa.utf8()),
            ("channel", pa.utf8()),
            ("intensity", pa.float32())
        ])
        metavalue_type = pa.struct([
            ("name", pa.utf8()),
            ("value", pa.utf8()),
            ("value_type", pa.utf8())
        ])

        data_dict = {}

        # QPX feature schema columns
        if should_include("sequence"):
            data_dict["sequence"] = pa.array(all_sequence, type=pa.utf8())
        if should_include("peptidoform"):
            data_dict["peptidoform"] = pa.array(all_peptidoform, type=pa.utf8())
        if should_include("modifications"):
            data_dict["modifications"] = pa.array(all_modifications, type=pa.list_(modification_type))
        if should_include("precursor_charge"):
            data_dict["precursor_charge"] = pa.array(all_precursor_charge, type=pa.int32())
        if should_include("calculated_mz"):
            data_dict["calculated_mz"] = pa.array(all_calculated_mz, type=pa.float32())
        if should_include("observed_mz"):
            data_dict["observed_mz"] = pa.array(all_observed_mz, type=pa.float32())
        if should_include("rt"):
            data_dict["rt"] = pa.array(all_rt, type=pa.float32())
        if should_include("posterior_error_probability"):
            data_dict["posterior_error_probability"] = pa.array(all_pep, type=pa.float64())
        if should_include("is_decoy"):
            data_dict["is_decoy"] = pa.array(all_is_decoy, type=pa.int32())
        if should_include("additional_scores"):
            data_dict["additional_scores"] = pa.array(all_additional_scores, type=pa.list_(score_type_struct))
        if should_include("predicted_rt"):
            data_dict["predicted_rt"] = pa.array(all_predicted_rt, type=pa.float32())
        if should_include("reference_file_name"):
            data_dict["reference_file_name"] = pa.array(all_reference_file, type=pa.utf8())
        if should_include("cv_params"):
            data_dict["cv_params"] = pa.array(all_cv_params, type=pa.utf8())
        if should_include("scan"):
            data_dict["scan"] = pa.array(all_scan, type=pa.utf8())
        if should_include("ion_mobility"):
            data_dict["ion_mobility"] = pa.array(all_ion_mobility, type=pa.float32())
        if should_include("start_ion_mobility"):
            data_dict["start_ion_mobility"] = pa.array(all_start_ion_mobility, type=pa.float32())
        if should_include("stop_ion_mobility"):
            data_dict["stop_ion_mobility"] = pa.array(all_stop_ion_mobility, type=pa.float32())
        if should_include("intensities"):
            data_dict["intensities"] = pa.array(all_intensities, type=pa.list_(intensity_type))
        if should_include("additional_intensities"):
            data_dict["additional_intensities"] = pa.array(all_additional_intensities, type=pa.utf8())
        if should_include("pg_accessions"):
            data_dict["pg_accessions"] = pa.array(all_pg_accessions, type=pa.list_(pa.utf8()))
        if should_include("anchor_protein"):
            data_dict["anchor_protein"] = pa.array(all_anchor_protein, type=pa.utf8())
        if should_include("unique"):
            data_dict["unique"] = pa.array(all_unique, type=pa.int32())
        if should_include("pg_global_qvalue"):
            data_dict["pg_global_qvalue"] = pa.array(all_pg_global_qvalue, type=pa.float64())
        if should_include("gg_accessions"):
            data_dict["gg_accessions"] = pa.array(all_gg_accessions, type=pa.list_(pa.utf8()))
        if should_include("gg_names"):
            data_dict["gg_names"] = pa.array(all_gg_names, type=pa.list_(pa.utf8()))
        if should_include("scan_reference_file_name"):
            data_dict["scan_reference_file_name"] = pa.array(all_scan_reference_file_name, type=pa.utf8())
        if should_include("rt_start"):
            data_dict["rt_start"] = pa.array(all_rt_start, type=pa.float32())
        if should_include("rt_stop"):
            data_dict["rt_stop"] = pa.array(all_rt_stop, type=pa.float32())

        # OpenMS-specific columns
        if should_include("quality"):
            data_dict["quality"] = pa.array(all_quality, type=pa.float32())
        if should_include("score"):
            data_dict["score"] = pa.array(all_score, type=pa.float64())
        if should_include("score_type"):
            data_dict["score_type"] = pa.array(all_score_type, type=pa.utf8())
        if should_include("spectrum_reference"):
            data_dict["spectrum_reference"] = pa.array(all_spectrum_reference, type=pa.utf8())
        if should_include("feature_metavalues"):
            all_feature_metavalues_str = [
                [{"name": mv["name"],
                  "value": str(mv["value"]) if mv["value"] is not None else None,
                  "value_type": mv["value_type"]}
                 for mv in mvs] for mvs in all_feature_metavalues
            ]
            data_dict["feature_metavalues"] = pa.array(all_feature_metavalues_str, type=pa.list_(metavalue_type))

        table = pa.Table.from_pydict(data_dict)

        if columns_set is not None:
            unknown = [c for c in columns if c not in data_dict]
            if unknown:
                warnings.warn(f"Unknown column(s) ignored in to_feature_arrow(): {unknown}. "
                              f"Use feature_columns() to see available columns.")

        # Sort by rt, observed_mz, precursor_charge
        sort_cols = ["rt", "observed_mz", "precursor_charge"]
        available_sort_cols = [(c, "ascending") for c in sort_cols if c in table.column_names]
        if available_sort_cols and table.num_rows > 0:
            import pyarrow.compute as pc
            indices = pc.sort_indices(table, sort_keys=available_sort_cols, null_placement="at_end")
            table = table.take(indices)

        return table

    def to_feature_df(self, **kwargs):
        """
        to_feature_df(self: ConsensusMap, **kwargs) -> pd.DataFrame

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export consensus features as DataFrame following QPX feature schema.

        Thin wrapper around to_feature_arrow().to_pandas().

        All keyword arguments are passed to to_feature_arrow().

        :return: DataFrame with QPX feature columns.
        :rtype: pd.DataFrame
        """
        return self.to_feature_arrow(**kwargs).to_pandas()

    def to_feature_parquet(self, path, compression='zstd', compression_level=None,
                            row_group_size=None, write_statistics=True, **kwargs):
        """
        to_feature_parquet(self: ConsensusMap, path: str, compression: str = 'zstd', compression_level: int = None, row_group_size: int = None, write_statistics: bool = True, **kwargs) -> None

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export consensus features to a Parquet file following QPX feature schema.

        :param path: Path to the output Parquet file.
        :type path: str
        :param compression: Compression algorithm ('none', 'snappy', 'gzip', 'lz4', 'zstd').
        :type compression: str
        :param compression_level: Compression level (depends on algorithm).
        :type compression_level: int
        :param row_group_size: Number of rows per row group.
        :type row_group_size: int
        :param write_statistics: Write column statistics. Default True.
        :type write_statistics: bool

        Additional kwargs are passed to to_feature_arrow().
        """
        try:
            import pyarrow.parquet as pq
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_feature_parquet(). "
                "Please install it with: pip install pyarrow"
            )

        table = self.to_feature_arrow(**kwargs)

        compression_map = {
            'none': None, 'snappy': 'snappy', 'gzip': 'gzip',
            'lz4': 'lz4', 'zstd': 'zstd',
        }
        compression_lower = compression.lower() if compression else 'none'
        if compression_lower not in compression_map:
            raise ValueError(f"Unsupported compression: {compression}. "
                           f"Use one of: {list(compression_map.keys())}")

        write_kwargs = {
            'compression': compression_map[compression_lower],
            'write_statistics': write_statistics,
        }
        if compression_level is not None and compression_map[compression_lower] in ('zstd', 'gzip'):
            write_kwargs['compression_level'] = compression_level
        if row_group_size is not None:
            write_kwargs['row_group_size'] = row_group_size

        with open(path, 'wb') as f:
            pq.write_table(table, f, **write_kwargs)

    def to_feature_qpx(self, qpx_version="1.0", creator="pyOpenMS", software_provider="OpenMS",
                        scan_format="scan", **kwargs):
        """
        to_feature_qpx(self: ConsensusMap, qpx_version: str = "1.0", creator: str = "pyOpenMS", software_provider: str = "OpenMS", scan_format: str = "scan", **kwargs) -> dict

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Export consensus features as QPX format structure with file metadata.

        :param qpx_version: QPX format version. Default "1.0".
        :type qpx_version: str
        :param creator: Creator name. Default "pyOpenMS".
        :type creator: str
        :param software_provider: Software provider. Default "OpenMS".
        :type software_provider: str
        :param scan_format: Controls scan column format.
        :type scan_format: str
        :return: Dict with 'file_metadata' and 'features' keys.
        :rtype: dict
        """
        import uuid as _uuid
        from datetime import datetime as _datetime

        df = self.to_feature_df(scan_format=scan_format, **kwargs)

        file_metadata = {
            "qpx_version": qpx_version,
            "creator": creator,
            "file_type": "feature",
            "creation_date": _datetime.now().isoformat(),
            "uuid": str(_uuid.uuid4()),
            "scan_format": scan_format,
            "software_provider": software_provider,
        }

        features = df.to_dict(orient="records")

        return {
            "file_metadata": file_metadata,
            "features": features
        }

    def feature_columns(self):
        """
        feature_columns(self: ConsensusMap) -> list

        **EXPERIMENTAL**: This method is experimental and subject to change.

        Return list of column names that to_feature_arrow() and to_feature_qpx() produce.

        :return: List of column name strings.
        :rtype: list
        """
        return [
            # QPX feature schema fields
            "sequence",
            "peptidoform",
            "modifications",
            "precursor_charge",
            "calculated_mz",
            "observed_mz",
            "rt",
            "posterior_error_probability",
            "is_decoy",
            "additional_scores",
            "predicted_rt",
            "reference_file_name",
            "cv_params",
            "scan",
            "ion_mobility",
            "start_ion_mobility",
            "stop_ion_mobility",
            "intensities",
            "additional_intensities",
            "pg_accessions",
            "anchor_protein",
            "unique",
            "pg_global_qvalue",
            "gg_accessions",
            "gg_names",
            "scan_reference_file_name",
            "rt_start",
            "rt_stop",
            # OpenMS-specific fields
            "quality",
            "score",
            "score_type",
            "spectrum_reference",
            "feature_metavalues",
        ]
