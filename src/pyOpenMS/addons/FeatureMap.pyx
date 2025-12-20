from UniqueIdInterface cimport setUniqueId as _setUniqueId
cimport numpy as np
import numpy as np




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

    def get_df_columns(self, columns='default', export_peptide_identifications=True):
        """
        get_df_columns(self: FeatureMap, columns: str = 'default', export_peptide_identifications: bool = True) -> List[str]

        Returns a list of column names that get_df() would produce.

        Useful for discovering available columns before export.

        Args:
            columns (str): 'default' for standard columns, 'all' to include all meta values.
            export_peptide_identifications (bool): Whether to include peptide ID columns.

        Returns:
            list: List of column name strings.

        Example:
            >>> fmap.get_df_columns()
            ['feature_id', 'peptide_sequence', 'charge', 'rt', 'mz', ...]
        """
        cols = ['feature_id']

        if export_peptide_identifications:
            cols.extend(['peptide_sequence', 'peptide_score', 'ID_filename', 'ID_native_id'])

        cols.extend(['charge', 'rt', 'mz', 'rt_start', 'rt_end', 'mz_start', 'mz_end', 'quality', 'intensity'])

        # Add meta values if 'all' requested
        if columns == 'all':
            meta_values = set()
            for f in self:
                mvs = []
                f.getKeys(mvs)
                for m in mvs:
                    meta_values.add(m.decode() if isinstance(m, bytes) else m)
            cols.extend(sorted(meta_values))

        return cols

    def _get_prot_id_filename_from_pep_id(self, pep_id):
        """
        _get_prot_id_filename_from_pep_id(self: FeatureMap, pep_id: PeptideIdentification) -> str

        Gets the primary MS run path of the ProteinIdentification linked with the given PeptideIdentification.

        Args:
        pep_id: PeptideIdentification

        Returns:
        str: primary MS run path (filename) of the ProteinIdentification with the same identifier as the given PeptideIdentification
        """
        for prot in self.getProteinIdentifications():
            if prot.getIdentifier() == pep_id.getIdentifier():
                filenames = []
                prot.getPrimaryMSRunPath(filenames)
                if filenames and filenames[0] != '':
                    return filenames[0]
        return 'unknown'

    def get_df(self, columns=None, meta_values=None, export_peptide_identifications=True):
        """
        get_df(self: FeatureMap, columns: Optional[List[str]] = None, meta_values: Optional[Union[List[str], str]] = None, export_peptide_identifications: bool = True) -> pd.DataFrame

        Generates a pandas DataFrame with information contained in the FeatureMap.

        Optionally the feature meta values and information for the assigned PeptideHit can be exported.

        :param columns: List of column names to include. If None,
                        includes all columns. Use get_df_columns() to discover available columns.
        :type columns: Optional[List[str]]

        :param meta_values: Meta values to include (None, [custom list of meta value names] or 'all')
        :type meta_values: Optional[Union[List[str], str]]

        :param export_peptide_identifications: Export sequence and score for best PeptideHit assigned to a feature.
                                               Additionally the ID_filename (file name of the corresponding ProteinIdentification)
                                               and the ID_native_id (spectrum ID of the corresponding Feature) are exported.
                                               They are also annotated as meta values when collecting all assigned PeptideIdentifications
                                               from a FeatureMap with FeatureMap.get_assigned_peptide_identifications().
                                               A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides)
                                               can be merged with the FeatureMap DataFrame with:
                                               merged_df = pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])
        :type export_peptide_identifications: bool

        :return: Feature information stored in a DataFrame
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed

        Example::

            # Get all columns
            df = fmap.get_df()

            # Discover available columns
            print(fmap.get_df_columns())

            # Get only specific columns
            df = fmap.get_df(columns=['feature_id', 'mz', 'rt', 'intensity'])
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for get_df(). "
                "Please install it with: pip install pandas"
            )
        # Common meta value types for numpy type mapping
        common_meta_value_types = {
            b'label': 'U50',
            b'spectrum_index': 'i',
            b'score_fit': 'f',
            b'score_correlation': 'f',
            b'FWHM': 'f',
            b'spectrum_native_id': 'U100',
            b'max_height': 'f',
            b'num_of_masstraces': 'i',
            b'masstrace_intensity': 'f',
            b'Group': 'U50',
            b'is_ungrouped_monoisotopic': 'i',
            b'leftWidth': 'f',
            b'rightWidth': 'f',
            b'total_xic': 'f',
            b'PeptideRef': 'U100',
            b'peak_apices_sum': 'f'
        }

        # Determine if peptide IDs are actually needed based on column selection
        pep_id_cols = {'peptide_sequence', 'peptide_score', 'ID_filename', 'ID_native_id'}
        if columns is not None:
            requested = set(columns)
            # Only export peptide IDs if explicitly requested or export_peptide_identifications is True
            # and at least one peptide column is requested
            need_pep_ids = export_peptide_identifications and len(requested & pep_id_cols) > 0
        else:
            need_pep_ids = export_peptide_identifications

        # get all possible meta value keys in a set
        if meta_values == 'all':
            meta_values = set()
            for f in self:
                mvs = []
                f.getKeys(mvs)
                for m in mvs:
                    meta_values.add(m)

        elif not meta_values: # if None, set to empty list
            meta_values = []

        def gen(fmap, fun):
            for f in fmap:
                yield from fun(f)

        def extract_meta_data(f):
            """Extracts feature meta data.

            Extracts information from a given feature with the requested meta values and, if requested,
            the sequence, score and ID_filename (primary MS run path of the linked ProteinIdentification)
            of the best PeptideHit (first) assigned to that feature.

            Args:
            f (Feature): feature from which to extract the meta data

            Yields:
            tuple: tuple containing feature information, peptide information (optional) and meta values (optional)
            """
            pep = f.getPeptideIdentifications()
            bb = f.getConvexHull().getBoundingBox2D()

            vals = [f.getMetaValue(m) if f.metaValueExists(m) else np.nan for m in meta_values]

            if need_pep_ids:
                if pep.size() > 0:
                    ID_filename = self._get_prot_id_filename_from_pep_id(pep[0])
                    # Check if spectrum_native_id meta value exists before accessing
                    spec_id = None
                    if f.metaValueExists('spectrum_native_id'):
                        spec_id = f.getMetaValue('spectrum_native_id')
                    hits = pep[0].getHits()
                    if len(hits) > 0:
                        besthit = hits[0]
                        pep_values = (besthit.getSequence().toString(), besthit.getScore(), ID_filename, spec_id)
                    else:
                        pep_values = (None, None, ID_filename, spec_id)
                else:
                    pep_values = (None, None, None, None)
            else:
                pep_values = ()

            yield tuple([f.getUniqueId()]) + pep_values + (f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals)

        cnt = self.size()

        mddtypes = [('feature_id', 'U100')]
        if need_pep_ids:
            mddtypes += [('peptide_sequence', 'U200'), ('peptide_score', 'f'), ('ID_filename', 'U100'), ('ID_native_id', 'U100')]
        mddtypes += [('charge', 'i4'), ('rt', np.dtype('double')), ('mz', np.dtype('double')), ('rt_start', np.dtype('double')), ('rt_end', np.dtype('double')),
                    ('mz_start', np.dtype('double')), ('mz_end', np.dtype('double')), ('quality', 'f'), ('intensity', 'f')]

        for meta_value in meta_values:
            if meta_value in common_meta_value_types:
                mddtypes.append((meta_value.decode(), common_meta_value_types[meta_value]))
            else:
                mddtypes.append((meta_value.decode(), 'U50'))

        mdarr = np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        df = pd.DataFrame(mdarr).set_index('feature_id')

        # Filter columns if requested
        if columns is not None:
            available_cols = [c for c in columns if c in df.columns or c == 'feature_id']
            # feature_id is the index, handle it specially
            if 'feature_id' not in available_cols:
                available_cols = [c for c in columns if c in df.columns]
            df = df[[c for c in available_cols if c in df.columns]]

        return df

    def to_arrow(self, columns=None, meta_values=None, export_peptide_identifications=True):
        """
        to_arrow(self: FeatureMap, columns: Optional[List[str]] = None, meta_values: Optional[Union[List[str], str]] = None, export_peptide_identifications: bool = True) -> pa.Table

        Returns an Apache Arrow Table with feature information.

        :param columns: List of column names to include. If None,
                        includes all columns. Use get_df_columns() to discover available columns.
        :type columns: Optional[List[str]]

        :param meta_values: Meta values to include (None, [custom list of meta value names] or 'all').
        :type meta_values: Optional[Union[List[str], str]]

        :param export_peptide_identifications: Export sequence and score for best PeptideHit
                                               assigned to a feature.
        :type export_peptide_identifications: bool

        :return: Arrow Table with feature data.
        :rtype: pyarrow.Table

        :raises ImportError: If pyarrow is not installed

        Example::

            # Get all columns
            table = fmap.to_arrow()

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
        df = self.get_df(columns=columns, meta_values=meta_values,
                        export_peptide_identifications=export_peptide_identifications)
        return pa.Table.from_pandas(df)

    def get_assigned_peptide_identifications(self):
        """
        get_assigned_peptide_identifications(self: FeatureMap) -> PeptideIdentificationList

        Generates a list with peptide identifications assigned to a feature.

        Adds 'ID_native_id' (feature spectrum id), 'ID_filename' (primary MS run path of corresponding ProteinIdentification)
        and 'feature_id' (unique ID of corresponding Feature) as meta values to the peptide hits.
        A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides) can be
        merged with the FeatureMap DataFrame with:
        merged_df = pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])

        Returns:
        PeptideIdentificationList: list of PeptideIdentification objects
        """
        from . import PeptideIdentificationList as _PeptideIdentificationList
        result = _PeptideIdentificationList()
        for f in self:
            for pep in f.getPeptideIdentifications():
                hits = []
                for hit in pep.getHits():
                    hit.setMetaValue('feature_id', str(f.getUniqueId()))
                    hit.setMetaValue('ID_filename', self._get_prot_id_filename_from_pep_id(pep))
                    if f.metaValueExists('spectrum_native_id'):
                        hit.setMetaValue('ID_native_id', f.getMetaValue('spectrum_native_id'))
                    else:
                        hit.setMetaValue('ID_native_id', 'unknown')
                    hits.append(hit)
                pep.setHits(hits)
                result.push_back(pep)
        return result
