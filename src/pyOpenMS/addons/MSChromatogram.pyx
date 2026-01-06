cimport numpy as np
import numpy as np




    def get_df_columns(self, columns='default', export_meta_values=True):
        """
        get_df_columns(self: MSChromatogram, columns: str = 'default', export_meta_values: bool = True) -> List[str]
        
        Returns a list of column names that get_df() would produce for this chromatogram.

        Useful for discovering available columns before export, especially when
        selecting specific columns for performance optimization.

        :param columns: 'default' for standard columns, 'all' for all available
                        columns including non-default ones (chromatogram_type, comment).
        :type columns: str
        :param export_meta_values: Whether to include meta value column names.
                                   Defaults to True.
        :type export_meta_values: bool
        :return: List of column name strings.
        :rtype: list

        Example::

            >>> # See default columns
            >>> cols = chrom.get_df_columns()
            ['rt', 'intensity', 'precursor_mz', ...]

            >>> # See ALL available columns
            >>> cols = chrom.get_df_columns('all')
            ['rt', 'intensity', ..., 'chromatogram_type', 'comment']

            >>> # Export everything
            >>> df = chrom.get_df(columns=chrom.get_df_columns('all'))
        """
        # Default columns (chromatogram_type and comment NOT included by default)
        cols = ['rt', 'intensity', 'precursor_mz', 'precursor_charge',
                'product_mz', 'native_id']

        # Add non-default columns if 'all' requested
        if columns == 'all':
            cols.extend(['chromatogram_type', 'comment'])

        if export_meta_values:
            mvs = []
            self.getKeys(mvs)
            for k in mvs:
                k_str = k.decode() if isinstance(k, bytes) else k
                cols.append(k_str)

        return cols

    def get_data_dict(self, columns=None, export_meta_values=True):
        """
        get_data_dict(self: MSChromatogram, columns: Optional[List[str]] = None, export_meta_values: bool = True) -> Dict[str, np.ndarray]
        
        Returns a dictionary of NumPy arrays with RT, intensities, and metadata.

        This method extracts chromatogram data including peaks, precursor/product info,
        and optional meta values into a dictionary format suitable for conversion to
        a pandas DataFrame.

        :param columns: List of column names to include. If None, includes
                        all default columns. Use get_df_columns('all') to see
                        all available columns.
        :type columns: Optional[List[str]]
        :param export_meta_values: Whether to include meta values in the output.
                                   Only applies when columns=None. Defaults to True.
        :type export_meta_values: bool
        :return: Dictionary with requested columns as keys and numpy arrays as values.
                 Default columns include:
                 - 'rt': numpy array of retention time values (float64)
                 - 'intensity': numpy array of intensity values (float32)
                 - 'precursor_mz': precursor m/z (float64)
                 - 'precursor_charge': precursor charge (uint16)
                 - 'product_mz': product m/z (float64)
                 - 'native_id': chromatogram native identifier
                 - Additional meta value columns (if export_meta_values=True)

                 Non-default columns (must be explicitly requested):
                 - 'chromatogram_type': type of chromatogram
                 - 'comment': chromatogram comment
        :rtype: dict

        Example::

            >>> # Get all columns (default)
            >>> data = chrom.get_data_dict()

            >>> # Get only specific columns for performance
            >>> data = chrom.get_data_dict(columns=['rt', 'intensity'])

            >>> # Get all available columns including non-defaults
            >>> all_cols = chrom.get_df_columns('all')
            >>> data = chrom.get_data_dict(columns=all_cols)
        """
        # Get peak data using existing optimized method
        cdef np.ndarray[np.float64_t, ndim=1] rts
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        rts, intensities = self.get_peaks()
        cnt = len(rts)

        # Determine which columns to include
        if columns is not None:
            requested = set(columns)
        else:
            requested = None  # None means include all defaults

        def want(col):
            """Check if a default column should be included."""
            return requested is None or col in requested

        def want_explicit(col):
            """Check if a non-default column is explicitly requested."""
            return requested is not None and col in requested

        data_dict = {}

        # Core peak data
        if want('rt'):
            data_dict['rt'] = rts
        if want('intensity'):
            data_dict['intensity'] = intensities

        # Precursor/Product info
        if want('precursor_mz'):
            data_dict['precursor_mz'] = np.full(cnt, self.getPrecursor().getMZ(), dtype=np.float64)
        if want('precursor_charge'):
            data_dict['precursor_charge'] = np.full(cnt, self.getPrecursor().getCharge(), dtype=np.uint16)
        if want('product_mz'):
            data_dict['product_mz'] = np.full(cnt, self.getProduct().getMZ(), dtype=np.float64)

        # Identifier
        if want('native_id'):
            data_dict['native_id'] = np.full(cnt, self.getNativeID(), dtype='object')

        # Non-default columns (only if explicitly requested)
        if want_explicit('chromatogram_type'):
            chrom_type = self.getChromatogramType()
            # Map enum value to name
            type_names = {
                0: 'MASS_CHROMATOGRAM',
                1: 'TOTAL_ION_CURRENT_CHROMATOGRAM',
                2: 'SELECTED_ION_CURRENT_CHROMATOGRAM',
                3: 'BASEPEAK_CHROMATOGRAM',
                4: 'SELECTED_ION_MONITORING_CHROMATOGRAM',
                5: 'SELECTED_REACTION_MONITORING_CHROMATOGRAM',
                6: 'ELECTROMAGNETIC_RADIATION_CHROMATOGRAM',
                7: 'ABSORPTION_CHROMATOGRAM',
                8: 'EMISSION_CHROMATOGRAM'
            }
            type_name = type_names.get(chrom_type, f'UNKNOWN_{chrom_type}')
            data_dict['chromatogram_type'] = np.full(cnt, type_name, dtype='object')

        if want_explicit('comment'):
            data_dict['comment'] = np.full(cnt, self.getComment(), dtype='object')

        # Meta values handling
        if requested is None and export_meta_values:
            mvs = []
            self.getKeys(mvs)
            for k in mvs:
                if not self.metaValueExists(k):
                    continue
                v = self.getMetaValue(k)
                k_str = k.decode() if isinstance(k, bytes) else k

                try:
                    # Check bool before int since bool is subclass of int in Python
                    if type(v) is type(True):
                        data_dict[k_str] = np.full(cnt, v, dtype=np.bool_)
                    elif isinstance(v, int):
                        data_dict[k_str] = np.full(cnt, v, dtype=np.int64)
                    elif isinstance(v, float):
                        data_dict[k_str] = np.full(cnt, v, dtype=np.float64)
                    elif isinstance(v, str):
                        data_dict[k_str] = np.full(cnt, v, dtype=f"U{max(len(v), 1)}")
                    else:
                        data_dict[k_str] = np.full(cnt, str(v), dtype='object')
                except Exception:
                    data_dict[k_str] = np.full(cnt, str(v), dtype='object')
        elif requested is not None:
            # Check if any requested columns are meta values
            mvs = []
            self.getKeys(mvs)
            mv_names = {(k.decode() if isinstance(k, bytes) else k): k for k in mvs}
            for col in requested:
                if col in mv_names:
                    k = mv_names[col]
                    if self.metaValueExists(k):
                        v = self.getMetaValue(k)
                        try:
                            if type(v) is type(True):
                                data_dict[col] = np.full(cnt, v, dtype=np.bool_)
                            elif isinstance(v, int):
                                data_dict[col] = np.full(cnt, v, dtype=np.int64)
                            elif isinstance(v, float):
                                data_dict[col] = np.full(cnt, v, dtype=np.float64)
                            elif isinstance(v, str):
                                data_dict[col] = np.full(cnt, v, dtype=f"U{max(len(v), 1)}")
                            else:
                                data_dict[col] = np.full(cnt, str(v), dtype='object')
                        except Exception:
                            data_dict[col] = np.full(cnt, str(v), dtype='object')

        return data_dict

    def get_df(self, columns=None, export_meta_values=True):
        """
        get_df(self: MSChromatogram, columns: Optional[List[str]] = None, export_meta_values: bool = True) -> pd.DataFrame

        Returns a pandas DataFrame representation of the MSChromatogram.

        This method converts the chromatogram data (peaks, metadata, precursor/product info)
        into a pandas DataFrame format.

        :param columns: List of column names to include. If None,
                        includes all default columns. Use get_df_columns()
                        to discover available columns.
        :type columns: Optional[List[str]]

        :param export_meta_values: Whether to include meta values. Only applies
                                   when columns=None. Defaults to True.
        :type export_meta_values: bool

        :return: DataFrame with requested columns. Default columns include:
                 rt, intensity, precursor_mz, precursor_charge, product_mz, native_id,
                 and additional meta value columns (if export_meta_values=True).
                 Non-default columns (must be explicitly requested): chromatogram_type, comment.
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed

        Example::

            # Get all default columns
            df = chrom.get_df()

            # Discover available columns
            print(chrom.get_df_columns())

            # Get only specific columns (faster)
            df = chrom.get_df(columns=['rt', 'intensity'])

            # Get all columns including non-defaults
            cols = chrom.get_df_columns('all')
            df = chrom.get_df(columns=cols)
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for get_df(). "
                "Please install it with: pip install pandas"
            )
        data_dict = self.get_data_dict(columns=columns, export_meta_values=export_meta_values)
        return pd.DataFrame(data_dict)

    def to_arrow(self, columns=None, export_meta_values=True):
        """
        to_arrow(self: MSChromatogram, columns: Optional[List[str]] = None, export_meta_values: bool = True) -> pa.Table

        Returns an Apache Arrow Table representation of the MSChromatogram.

        This method converts the chromatogram data (peaks, metadata, precursor/product info)
        into an Arrow Table format for efficient data interchange.

        :param columns: List of column names to include. If None,
                        includes all default columns. Use get_df_columns()
                        to discover available columns.
        :type columns: Optional[List[str]]

        :param export_meta_values: Whether to include meta values. Only applies
                                   when columns=None. Defaults to True.
        :type export_meta_values: bool

        :return: Arrow Table with requested columns.
        :rtype: pyarrow.Table

        :raises ImportError: If pyarrow is not installed

        Example::

            # Get all default columns
            table = chrom.to_arrow()

            # Get only specific columns (faster)
            table = chrom.to_arrow(columns=['rt', 'intensity'])

            # Convert to pandas (zero-copy with pandas 2.0+)
            df = table.to_pandas()
        """
        try:
            import pyarrow as pa
        except ImportError:
            raise ImportError(
                "pyarrow is required for to_arrow(). "
                "Please install it with: pip install pyarrow"
            )
        data_dict = self.get_data_dict(columns=columns, export_meta_values=export_meta_values)
        return pa.Table.from_pydict(data_dict)

    def get_peaks(self):
        """
        get_peaks(self: MSChromatogram) -> Tuple[np.ndarray, np.ndarray]
        
        Get RT and intensity arrays as numpy arrays.
        """
        cdef _MSChromatogram * chrom_ = self.inst.get()

        cdef unsigned int n = chrom_.size()
        cdef np.ndarray[np.float64_t, ndim=1] rts
        rts = np.zeros( (n,), dtype=np.float64)
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.zeros( (n,), dtype=np.float32)
        cdef _ChromatogramPeak p

        cdef libcpp_vector[_ChromatogramPeak].iterator it = chrom_.begin()
        cdef int i = 0
        while it != chrom_.end():
            rts[i] = deref(it).getRT()
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return rts, intensities

    def set_peaks(self, peaks):
        """
        set_peaks(self: MSChromatogram, peaks: Union[Tuple[np.ndarray, np.ndarray], Tuple[List, List]]) -> None
        
        Set RT and intensity data from arrays or lists.
        """
        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (rt and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (rt and intensity vector)"

        rts, intensities = peaks
        assert len(rts) == len(intensities), "Input vectors for set_peaks need to have the same length (rt and intensity vector)"

        # Select which function to use for set_peaks:
        # If we have numpy arrays, it helps to use optimized functions
        if isinstance(rts, np.ndarray) and isinstance(intensities, np.ndarray) and \
          rts.dtype == np.float64 and intensities.dtype == np.float32 and \
          rts.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_df(rts, intensities)
        elif isinstance(rts, np.ndarray) and isinstance(intensities, np.ndarray) and \
          rts.dtype == np.float64 and intensities.dtype == np.float64 and \
          rts.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_dd(rts, intensities)
        else:
            self._set_peaks_orig(rts, intensities)



    def _set_peaks_fast_dd(self, np.ndarray[double, ndim=1, mode="c"] data_rt not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):
        """
        _set_peaks_fast_dd(self: MSChromatogram, data_rt: np.ndarray, data_i: np.ndarray) -> None
        
        Internal method to set peaks from double/double numpy arrays.
        """
        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.resize(0) # empty vector, keep meta data and data arrays
        chrom_.reserve(<int>len(data_rt)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef double intensity
        cdef int N
        N = len(data_rt)

        for i in range(N):
            rt = data_rt[i]
            intensity = data_i[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()


    def _set_peaks_fast_df(self, np.ndarray[double, ndim=1, mode="c"] data_rt not None, np.ndarray[float, ndim=1, mode="c"] data_i not None):
        """
        _set_peaks_fast_df(self: MSChromatogram, data_rt: np.ndarray, data_i: np.ndarray) -> None
        
        Internal method to set peaks from double/float numpy arrays.
        """
        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.resize(0) # empty vector, keep meta data and data arrays
        chrom_.reserve(<int>len(data_rt)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef float intensity
        cdef int N
        N = len(data_rt)

        for i in range(N):
            rt = data_rt[i]
            intensity = data_i[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()


    def _set_peaks_orig(self, rts, intensities):
        """
        _set_peaks_orig(self: MSChromatogram, rts: Union[List, np.ndarray], intensities: Union[List, np.ndarray]) -> None
        
        Internal method to set peaks from generic lists or arrays.
        """
        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.resize(0) # empty vector, keep meta data and data arrays
        chrom_.reserve(<int>len(rts)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef float intensity
        cdef int N
        N = len(rts)

        for i in range(N):
            rt = rts[i]
            intensity = intensities[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()

    def __len__(self):
        """
        __len__(self: MSChromatogram) -> int
        
        Return the number of peaks in the chromatogram.
        """
        return self.inst.get().size()

    def __str__(self):
        """
        __str__(self: MSChromatogram) -> str
        
        Return a string representation of the MSChromatogram object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: MSChromatogram) -> str
        
        Return a string representation of the MSChromatogram object.

        Returns key properties in a readable format:
        MSChromatogram(name='transition_1', num_peaks=100, rt_range=[0.0, 3600.0])
        """
        cdef size_t num_peaks = self.inst.get().size()

        parts = []

        # Add name if set
        name = self.getName()
        if name:
            name_str = name.decode('utf-8') if isinstance(name, bytes) else str(name)
            if name_str:
                parts.append(f"name='{name_str}'")

        parts.append(f"num_peaks={num_peaks}")

        # Add RT range if there are peaks
        if num_peaks > 0:
            rts, _ = self.get_peaks()
            parts.append(f"rt_range=[{rts[0]:.2f}, {rts[-1]:.2f}]")

        return f"MSChromatogram({', '.join(parts)})"
