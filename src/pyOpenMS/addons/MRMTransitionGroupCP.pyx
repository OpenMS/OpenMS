cimport numpy as np
import numpy as np




    def get_chromatogram_df_columns(self, columns='default', export_meta_values=True):
        """
        get_chromatogram_df_columns(self: MRMTransitionGroupCP, columns: str = 'default', export_meta_values: bool = True) -> List[str]

        Returns a list of column names that get_chromatogram_df() would produce.

        Args:
            columns (str): 'default' for standard columns, 'all' for all available columns.
            export_meta_values (bool): Whether to include meta value columns.

        Returns:
            list: List of column name strings.
        """
        # Use the first chromatogram to get columns (all should have same structure)
        chroms = self.getChromatograms()
        if chroms:
            return chroms[0].get_df_columns(columns=columns, export_meta_values=export_meta_values)
        return ['rt', 'intensity', 'precursor_mz', 'precursor_charge', 'product_mz', 'native_id']

    def get_feature_df_columns(self, columns='default'):
        """
        get_feature_df_columns(self: MRMTransitionGroupCP, columns: str = 'default') -> List[str]

        Returns a list of column names that get_feature_df() would produce.

        Args:
            columns (str): 'default' for core columns, 'all' to include all meta values.

        Returns:
            list: List of column name strings.
        """
        cols = ['feature_id', 'rt', 'intensity', 'quality']

        if columns == 'all':
            meta_values = set()
            for f in self.getFeatures():
                mvs = []
                f.getKeys(mvs)
                for m in mvs:
                    meta_values.add(m.decode() if isinstance(m, bytes) else m)
            cols.extend(sorted(meta_values))

        return cols

    def get_chromatogram_df(self, columns=None, export_meta_values=True):
        """
        get_chromatogram_df(self: MRMTransitionGroupCP, columns: Optional[List[str]] = None, export_meta_values: bool = True) -> pd.DataFrame

        Returns a DataFrame representation of the Chromatograms stored in MRMTransitionGroupCP.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all default columns.
            export_meta_values (bool): Whether to export meta values. Only applies
                                       when columns=None.

        Returns:
            pd.DataFrame: DataFrame representation of the chromatograms.

        Example:
            >>> # Get all default columns
            >>> df = mrm.get_chromatogram_df()

            >>> # Discover available columns
            >>> print(mrm.get_chromatogram_df_columns())

            >>> # Get only specific columns
            >>> df = mrm.get_chromatogram_df(columns=['rt', 'intensity'])
        """
        import pandas as pd
        chroms = self.getChromatograms()
        out = [c.get_df(columns=columns, export_meta_values=export_meta_values) for c in chroms]
        if out:
            return pd.concat(out, ignore_index=True)
        return pd.DataFrame()

    def get_feature_df(self, columns=None, meta_values=None):
        """
        get_feature_df(self: MRMTransitionGroupCP, columns: Optional[List[str]] = None, meta_values: Optional[Union[List[str], str]] = None) -> pd.DataFrame

        Returns a DataFrame representation of the Features stored in MRMTransitionGroupCP.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all columns. Use get_feature_df_columns()
                                   to discover available columns.
            meta_values: meta values to include (None, [custom list of meta value names] or 'all')

        Returns:
            pd.DataFrame: DataFrame representation of the Features.

        Example:
            >>> # Get all columns
            >>> df = mrm.get_feature_df()

            >>> # Discover available columns
            >>> print(mrm.get_feature_df_columns())

            >>> # Get only specific columns
            >>> df = mrm.get_feature_df(columns=['feature_id', 'rt', 'intensity'])
        """
        import pandas as pd
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

        def gen(features, fun):
            for f in features:
                yield from fun(f)

        def extract_meta_data(f):
            """Extracts feature meta data.

            Extracts information from a given feature with the requested meta values.

            Parameters:
            f (MRMFeature): feature from which to extract the meta data

            Yields:
            tuple: tuple containing feature information and meta values (optional)
            """
            vals = [f.getMetaValue(m) if f.metaValueExists(m) else np.nan for m in meta_values_list]

            yield tuple((f.getUniqueId(), f.getRT(), f.getIntensity(), f.getOverallQuality(), *vals))

        # get all possible meta value keys in a set
        features = self.getFeatures()
        mddtypes = [('feature_id', np.dtype('uint64')), ('rt', 'f'), ('intensity', 'f'), ('quality', 'f')]

        if meta_values is not None:
            # Add all possible meta values to meta_value array if 'all' is passed
            if meta_values == 'all':
                meta_values_list = set()
                for f in features:
                    mvs = []
                    f.getKeys(mvs)
                    for m in mvs:
                        meta_values_list.add(m)
                meta_values_list = list(meta_values_list)
            else:
                meta_values_list = list(meta_values)

            # Add meta_values to mddtypes
            for meta_value in meta_values_list:
                if meta_value in common_meta_value_types:
                    mddtypes.append((meta_value.decode() if isinstance(meta_value, bytes) else meta_value, common_meta_value_types.get(meta_value, 'U50')))
                else:
                    mddtypes.append((meta_value.decode() if isinstance(meta_value, bytes) else meta_value, 'U50'))
        else:
            meta_values_list = []

        mdarr = np.fromiter(iter=gen(features, extract_meta_data), dtype=mddtypes, count=len(features))

        df = pd.DataFrame(mdarr).set_index('feature_id')

        # Filter columns if requested
        if columns is not None:
            available_cols = [c for c in columns if c in df.columns or c == 'feature_id']
            if 'feature_id' not in available_cols:
                available_cols = [c for c in columns if c in df.columns]
            df = df[[c for c in available_cols if c in df.columns]]

        return df

    def chromatograms_to_arrow(self, columns=None, export_meta_values=True):
        """
        chromatograms_to_arrow(self: MRMTransitionGroupCP, columns: Optional[List[str]] = None, export_meta_values: bool = True) -> pa.Table

        Returns an Apache Arrow Table representation of the Chromatograms.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all default columns.
            export_meta_values (bool): Whether to export meta values.

        Returns:
            pyarrow.Table: Arrow Table with chromatogram data.

        Example:
            >>> table = mrm.chromatograms_to_arrow()
            >>> df = table.to_pandas()
        """
        import pyarrow as pa
        df = self.get_chromatogram_df(columns=columns, export_meta_values=export_meta_values)
        return pa.Table.from_pandas(df)

    def features_to_arrow(self, columns=None, meta_values=None):
        """
        features_to_arrow(self: MRMTransitionGroupCP, columns: Optional[List[str]] = None, meta_values: Optional[Union[List[str], str]] = None) -> pa.Table

        Returns an Apache Arrow Table representation of the Features.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all columns.
            meta_values: Meta values to include (None, [custom list] or 'all').

        Returns:
            pyarrow.Table: Arrow Table with feature data.

        Example:
            >>> table = mrm.features_to_arrow()
            >>> df = table.to_pandas()
        """
        import pyarrow as pa
        df = self.get_feature_df(columns=columns, meta_values=meta_values)
        return pa.Table.from_pandas(df)
