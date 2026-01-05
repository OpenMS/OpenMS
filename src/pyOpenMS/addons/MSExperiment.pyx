cimport numpy as np
import numpy as np


    def get2DPeakData(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty] get2DPeakData(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""

        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[libcpp_vector[float]] mz
        cdef libcpp_vector[libcpp_vector[float]] inty
        exp_.get2DPeakDataPerSpectrum(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty)

        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)

        cdef np.ndarray all_mz = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_inty = np.empty(rt.size(), dtype=object)
        cdef ArrayWrapperFloat mz_wrap
        cdef ArrayWrapperFloat inty_wrap

        cdef unsigned int i
        for i in range(0, mz.size()):
            mz_wrap = ArrayWrapperFloat()
            inty_wrap = ArrayWrapperFloat()
            mz_wrap.set_data(mz[i])
            inty_wrap.set_data(inty[i])
            all_mz[i] = np.frombuffer(mz_wrap)
            all_inty[i] = np.frombuffer(inty_wrap)

        return (np.frombuffer(rt_wrap), all_mz, all_inty)



    def get2DPeakDataLong(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty] get2DPeakDataLong(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[float] mz
        cdef libcpp_vector[float] inty
        exp_.get2DPeakData(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty)
       
        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat mz_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat inty_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)
        mz_wrap.set_data(mz)
        inty_wrap.set_data(inty)

        return (np.asarray(rt_wrap), np.asarray(mz_wrap), np.asarray(inty_wrap))
    
    def get2DPeakDataIM(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty, np.array[float] ion_mobility] get2DPeakDataIM(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[libcpp_vector[float]] mz
        cdef libcpp_vector[libcpp_vector[float]] inty
        cdef libcpp_vector[libcpp_vector[float]] ion_mobility
        exp_.get2DPeakDataIMPerSpectrum(min_rt, max_rt, min_mz, max_mz,  ms_level, rt, mz, inty, ion_mobility)

        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)

        cdef np.ndarray all_mz = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_inty = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_ion = np.empty(rt.size(), dtype=object)
        cdef ArrayWrapperFloat mz_wrap
        cdef ArrayWrapperFloat inty_wrap
        cdef ArrayWrapperFloat ion_mobility_wrap

        cdef unsigned int i
        for i in range(0, mz.size()):
            mz_wrap = ArrayWrapperFloat()
            inty_wrap = ArrayWrapperFloat()
            ion_mobility_wrap = ArrayWrapperFloat()
            mz_wrap.set_data(mz[i])
            inty_wrap.set_data(inty[i])
            ion_mobility_wrap.set_data(ion_mobility[i])
            all_mz[i] = np.frombuffer(mz_wrap)
            all_inty[i] = np.frombuffer(inty_wrap)
            all_ion[i] = np.frombuffer(ion_mobility_wrap)

        return (np.frombuffer(rt_wrap), all_mz, all_inty, all_ion)

    def get2DPeakDataIMLong(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty, np.array[float] ion_mobility] get2DPeakDataIMLong(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[float] mz
        cdef libcpp_vector[float] inty
        cdef libcpp_vector[float] ion_mobility
        exp_.get2DPeakDataIM(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty, ion_mobility)
       
        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat mz_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat inty_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat ion_mobility_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)
        mz_wrap.set_data(mz)
        inty_wrap.set_data(inty)
        ion_mobility_wrap.set_data(ion_mobility)

        return (np.asarray(rt_wrap), np.asarray(mz_wrap), np.asarray(inty_wrap), np.asarray(ion_mobility_wrap))

    def getMSLevels(self):
        """Cython signature: list[int] getMSLevels()"""
        cdef libcpp_vector[unsigned int] _r = self.inst.get().getMSLevels()
        cdef libcpp_vector[unsigned int].iterator it__r = _r.begin()
        cdef list result = []
        while it__r != _r.end():
            result.append(deref(it__r))
            inc(it__r)
        return result

    def getChromatogram(self, id_):
        """Cython signature: `MSChromatogram getChromatogram(size_t id_)`"""
        assert isinstance(id_, int), 'arg id_ wrong type'
        assert id_ < self.getNrChromatograms(), 'Requested chromatogram %s does not exist, there are only %s chromatograms' % (id_, self.getNrChromatograms() )
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogram((<size_t>id_)))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result

    def getSpectrum(self, id_):
        """Cython signature: `MSSpectrum getSpectrum(size_t id_)`"""
        assert isinstance(id_, int), 'arg id_ wrong type'
        assert id_ < self.getNrSpectra(), 'Requested spectrum %s does not exist, there are only %s spectra' % (id_, self.getNrSpectra() )
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrum((<size_t>id_)))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result

    def __len__(self):
        """Return the number of spectra in the experiment."""
        return self.inst.get().size()

    def __str__(self):
        """
        __str__(self: MSExperiment) -> str
        
        Return a string representation of the MSExperiment object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: MSExperiment) -> str
        
        Return a string representation of the MSExperiment object.

        Returns key properties in a readable format:
        MSExperiment(num_spectra=1000, num_chromatograms=10, ms_levels=[1, 2], rt_range=[0.0, 3600.0])
        """
        cdef size_t num_spectra = self.getNrSpectra()
        cdef size_t num_chromatograms = self.getNrChromatograms()

        parts = []
        parts.append(f"num_spectra={num_spectra}")
        parts.append(f"num_chromatograms={num_chromatograms}")

        # Add MS levels
        ms_levels = self.getMSLevels()
        if ms_levels:
            parts.append(f"ms_levels={ms_levels}")

        # Add RT range if there are spectra
        if num_spectra > 0:
            min_rt = self.getMinRT()
            max_rt = self.getMaxRT()
            parts.append(f"rt_range=[{min_rt:.2f}, {max_rt:.2f}]")

        return f"MSExperiment({', '.join(parts)})"

    def get_df_columns(self, long_format=False):
        """
        get_df_columns(self: MSExperiment, long_format: bool = False) -> List[str]

        Returns a list of column names that get_df() would produce.

        Useful for discovering available columns before export.

        :param long_format: If True, returns columns for long format.
                            If False, returns columns for compact format.
        :type long_format: bool
        :return: List of column name strings.
        :rtype: list

        Example::

            >>> exp.get_df_columns(long_format=True)
            ['rt', 'mz', 'intensity', 'ms_level']

            >>> exp.get_df_columns(long_format=False)
            ['rt', 'ms_level', 'mz_array', 'intensity_array']
        """
        if long_format:
            return ['rt', 'mz', 'intensity', 'ms_level']
        else:
            return ['rt', 'ms_level', 'mz_array', 'intensity_array']

    def get_arrow_columns(self, data='spectra', format='long',
                          include_precursor_info=True, include_ion_mobility=True):
        """
        get_arrow_columns(self: MSExperiment, data: str = 'spectra', format: str = 'long',
                          include_precursor_info: bool = True, include_ion_mobility: bool = True) -> List[str]

        Returns a list of column names that to_arrow() would produce with the given parameters.

        Useful for discovering available columns before export, especially for column selection.

        :param data: Type of data to export: 'spectra', 'chromatograms', or 'both'.
        :type data: str
        :param format: Output format: 'long' (one row per peak) or 'semi_wide' (one row per spectrum/chromatogram).
        :type format: str
        :param include_precursor_info: Include precursor columns (precursor_mz, precursor_charge, etc.).
        :type include_precursor_info: bool
        :param include_ion_mobility: Include ion mobility column (spectra only).
        :type include_ion_mobility: bool
        :return: List of column name strings, or dict of lists if data='both'.
        :rtype: Union[List[str], Dict[str, List[str]]]

        Example::

            >>> exp.get_arrow_columns(data='spectra', format='long')
            ['mz', 'intensity', 'rt', 'spectrum_index', 'ms_level', 'native_id', ...]

            >>> exp.get_arrow_columns(data='chromatograms', format='semi_wide')
            ['rt', 'intensity', 'chromatogram_index', 'native_id', ...]

            >>> exp.get_arrow_columns(data='both')
            {'spectra': [...], 'chromatograms': [...]}
        """
        # Column order matches C++ ArrowExport implementation
        spectra_long_cols = ['mz', 'intensity', 'rt']
        spectra_semi_wide_cols = ['mz', 'intensity', 'rt']

        # ion_mobility comes after rt in C++ (before spectrum_index)
        if include_ion_mobility:
            spectra_long_cols.append('ion_mobility')
            spectra_semi_wide_cols.append('ion_mobility')

        spectra_long_cols.extend(['spectrum_index', 'ms_level', 'native_id'])
        spectra_semi_wide_cols.extend(['spectrum_index', 'ms_level', 'native_id'])

        if include_precursor_info:
            precursor_cols = ['precursor_mz', 'precursor_charge', 'precursor_intensity',
                              'isolation_lower', 'isolation_upper']
            spectra_long_cols.extend(precursor_cols)
            spectra_semi_wide_cols.extend(precursor_cols)

        chrom_long_cols = ['rt', 'intensity', 'chromatogram_index', 'native_id',
                          'precursor_mz', 'product_mz']
        chrom_semi_wide_cols = ['rt', 'intensity', 'chromatogram_index', 'native_id',
                               'precursor_mz', 'product_mz']

        if data == 'spectra':
            return spectra_long_cols if format == 'long' else spectra_semi_wide_cols
        elif data == 'chromatograms':
            return chrom_long_cols if format == 'long' else chrom_semi_wide_cols
        elif data == 'both':
            return {
                'spectra': spectra_long_cols if format == 'long' else spectra_semi_wide_cols,
                'chromatograms': chrom_long_cols if format == 'long' else chrom_semi_wide_cols
            }
        else:
            raise ValueError(f"data must be 'spectra', 'chromatograms', or 'both', got '{data}'")

    def get_df(self, columns=None, ms_levels=None, long_format=False):
        """
        get_df(self: MSExperiment, columns: Optional[List[str]] = None, ms_levels: List[int] = [], long_format: bool = False) -> pd.DataFrame

        Generates a pandas DataFrame with all peaks in the MSExperiment

        :param columns: List of column names to include. If None,
                        includes all columns. Use get_df_columns() to discover available columns.
        :type columns: Optional[List[str]]

        :param ms_levels: Get only spectra with the given MS levels. Default is an empty list,
                          which means all MS levels will be included.
        :type ms_levels: List[int]

        :param long_format: Set to True if you want to have a long/expanded/melted dataframe
                            with one row per peak. Faster but replicated RT information.
                            If False, returns rows in the style: rt, np.array(mz), np.array(int)
        :type long_format: bool

        :return: Feature information stored in a DataFrame
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed

        Example::

            # Get all columns
            df = exp.get_df()

            # Discover available columns
            print(exp.get_df_columns())

            # Get only specific columns
            df = exp.get_df(columns=['rt', 'mz', 'intensity'], long_format=True)
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for get_df(). "
                "Please install it with: pip install pandas"
            )
        if ms_levels is None:
            ms_levels = []
        self.updateRanges()
        if not ms_levels:
            ms_levels = self.getMSLevels()
        if long_format:
            cols = ["rt", "mz", "intensity"]
            dfs = []
            for ms_level in ms_levels:
                spectraarrs2d = self.get2DPeakDataLong(self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ(), ms_level)
                df = pd.DataFrame(dict(zip(cols, spectraarrs2d)))
                df["ms_level"] = ms_level
                dfs.append(df)
            df = pd.concat(dfs, ignore_index=True)
        else:
            cols = ["rt", "ms_level", "mz_array", "intensity_array"]
            df = pd.DataFrame(data=((spec.getRT(), spec.getMSLevel(), *spec.get_peaks()) for spec in self if spec.getMSLevel() in ms_levels), columns=cols)

        # Filter columns if requested
        if columns is not None:
            available_cols = [c for c in columns if c in df.columns]
            df = df[available_cols]

        return df

    def to_arrow(self, data='spectra', format='long', columns=None, ms_levels=None,
                  min_rt=None, max_rt=None, min_mz=None, max_mz=None,
                  include_precursor_info=True, include_ion_mobility=True,
                  long_format=None):
        """
        to_arrow(self: MSExperiment, data: str = 'spectra', format: str = 'long', ...) -> Union[pa.Table, Dict[str, pa.Table]]

        Returns an Apache Arrow Table with spectra and/or chromatogram data from the MSExperiment.

        This is a unified interface that supports exporting both spectra and chromatograms
        in either long format (one row per peak/point) or semi-wide format (one row per
        spectrum/chromatogram with list columns for m/z and intensity arrays).

        :param data: Type of data to export:
                     - 'spectra': Export spectrum peak data (default)
                     - 'chromatograms': Export chromatogram data
                     - 'both': Export both as a dict {'spectra': Table, 'chromatograms': Table}
        :type data: str

        :param format: Output format:
                       - 'long': One row per peak/point (default). Scalars are repeated per peak.
                       - 'semi_wide': One row per spectrum/chromatogram. m/z and intensity are list arrays.
        :type format: str

        :param columns: List of column names to include. If None, includes all available columns.
                        Use get_arrow_columns() to discover available columns.
        :type columns: Optional[List[str]]

        :param ms_levels: Filter spectra by MS levels. Default None includes all levels.
                          Only applies when data='spectra' or 'both'.
        :type ms_levels: Optional[List[int]]

        :param min_rt: Minimum retention time filter. Default None means no filter.
        :type min_rt: Optional[float]

        :param max_rt: Maximum retention time filter. Default None means no filter.
        :type max_rt: Optional[float]

        :param min_mz: Minimum m/z filter (spectra only). Default None means no filter.
        :type min_mz: Optional[float]

        :param max_mz: Maximum m/z filter (spectra only). Default None means no filter.
        :type max_mz: Optional[float]

        :param include_precursor_info: Include precursor columns (precursor_mz, precursor_charge,
                                       precursor_intensity, isolation_lower, isolation_upper).
                                       Default True.
        :type include_precursor_info: bool

        :param include_ion_mobility: Include ion_mobility column if data is present. Default True.
                                     Only applies when data='spectra' or 'both'.
        :type include_ion_mobility: bool

        :param long_format: DEPRECATED. Use format='long' or format='semi_wide' instead.
                            If provided, overrides format parameter for backward compatibility.
        :type long_format: Optional[bool]

        :return: Arrow Table with the requested data.
                 If data='both', returns dict {'spectra': Table, 'chromatograms': Table}.
        :rtype: Union[pyarrow.Table, Dict[str, pyarrow.Table]]

        :raises ImportError: If pyarrow is not installed
        :raises ValueError: If data parameter is invalid

        Example::

            # Export spectra in long format (default)
            table = exp.to_arrow()

            # Export chromatograms in semi-wide format
            chrom_table = exp.to_arrow(data='chromatograms', format='semi_wide')

            # Export both spectra and chromatograms
            tables = exp.to_arrow(data='both')
            spectra_table = tables['spectra']
            chrom_table = tables['chromatograms']

            # Filter by MS level and RT range
            table = exp.to_arrow(ms_levels=[2], min_rt=100.0, max_rt=500.0)

            # Select specific columns
            table = exp.to_arrow(columns=['mz', 'intensity', 'rt', 'precursor_mz'])

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

        # Handle deprecated long_format parameter for backward compatibility
        if long_format is not None:
            import warnings
            warnings.warn(
                "long_format parameter is deprecated. Use format='long' or format='semi_wide' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            format = 'long' if long_format else 'semi_wide'

        # Validate parameters
        if data not in ('spectra', 'chromatograms', 'both'):
            raise ValueError(f"data must be 'spectra', 'chromatograms', or 'both', got '{data}'")
        if format not in ('long', 'semi_wide'):
            raise ValueError(f"format must be 'long' or 'semi_wide', got '{format}'")

        self.updateRanges()

        # Set default RT/MZ ranges (0 means no filter in C++ API)
        _min_rt = 0.0 if min_rt is None else min_rt
        _max_rt = 0.0 if max_rt is None else max_rt
        _min_mz = 0.0 if min_mz is None else min_mz
        _max_mz = 0.0 if max_mz is None else max_mz

        # Try zero-copy C++ path first (only available when built with WITH_PARQUET)
        try:
            from pyopenms._arrow_zerocopy import spectra_to_arrow, chromatograms_to_arrow
            _use_zerocopy = True
        except ImportError:
            _use_zerocopy = False

        if _use_zerocopy:
            # Use zero-copy C++ implementation
            result = {}

            if data in ('spectra', 'both'):
                result['spectra'] = spectra_to_arrow(
                    self,
                    format=format,
                    ms_levels=ms_levels,
                    min_rt=_min_rt, max_rt=_max_rt,
                    min_mz=_min_mz, max_mz=_max_mz,
                    columns=columns,
                    include_precursor_info=include_precursor_info,
                    include_ion_mobility=include_ion_mobility
                )

            if data in ('chromatograms', 'both'):
                result['chromatograms'] = chromatograms_to_arrow(
                    self,
                    format=format,
                    min_rt=_min_rt, max_rt=_max_rt,
                    columns=columns
                )

            if data == 'both':
                return result
            elif data == 'spectra':
                return result['spectra']
            else:
                return result['chromatograms']

        # Fallback to Python implementation
        # Handle empty experiment case
        num_spectra = self.getNrSpectra()
        num_chromatograms = self.getNrChromatograms()

        if num_spectra == 0 and num_chromatograms == 0:
            # Return empty tables for empty experiment with full schema
            # Column order matches C++ ArrowExport implementation
            result = {}
            if data in ('spectra', 'both'):
                spectra_dict = {
                    'mz': pa.array([], type=pa.float64()),
                    'intensity': pa.array([], type=pa.float32()),
                    'rt': pa.array([], type=pa.float32()),
                }
                # ion_mobility comes after rt (before spectrum_index) to match C++
                if include_ion_mobility:
                    spectra_dict['ion_mobility'] = pa.array([], type=pa.float32())
                spectra_dict['spectrum_index'] = pa.array([], type=pa.uint32())
                spectra_dict['ms_level'] = pa.array([], type=pa.uint8())
                spectra_dict['native_id'] = pa.array([], type=pa.utf8())
                if include_precursor_info:
                    spectra_dict['precursor_mz'] = pa.array([], type=pa.float64())
                    spectra_dict['precursor_charge'] = pa.array([], type=pa.int16())
                    spectra_dict['precursor_intensity'] = pa.array([], type=pa.float32())
                    spectra_dict['isolation_lower'] = pa.array([], type=pa.float64())
                    spectra_dict['isolation_upper'] = pa.array([], type=pa.float64())
                # Filter columns if requested
                if columns is not None:
                    spectra_dict = {k: v for k, v in spectra_dict.items() if k in columns}
                result['spectra'] = pa.Table.from_pydict(spectra_dict)

            if data in ('chromatograms', 'both'):
                chrom_dict = {
                    'rt': pa.array([], type=pa.float64()),
                    'intensity': pa.array([], type=pa.float32()),
                    'chromatogram_index': pa.array([], type=pa.uint32()),
                    'native_id': pa.array([], type=pa.utf8()),
                    'precursor_mz': pa.array([], type=pa.float64()),
                    'product_mz': pa.array([], type=pa.float64()),
                }
                # Filter columns if requested
                if columns is not None:
                    chrom_dict = {k: v for k, v in chrom_dict.items() if k in columns}
                result['chromatograms'] = pa.Table.from_pydict(chrom_dict)

            if data == 'both':
                return result
            elif data == 'spectra':
                return result['spectra']
            else:
                return result['chromatograms']

        # Set default RT/MZ ranges for Python path
        if min_rt is None and num_spectra > 0:
            min_rt = self.getMinRT()
        elif min_rt is None:
            min_rt = 0.0
        if max_rt is None and num_spectra > 0:
            max_rt = self.getMaxRT()
        elif max_rt is None:
            max_rt = float('inf')
        if min_mz is None and num_spectra > 0:
            min_mz = self.getMinMZ()
        elif min_mz is None:
            min_mz = 0.0
        if max_mz is None and num_spectra > 0:
            max_mz = self.getMaxMZ()
        elif max_mz is None:
            max_mz = float('inf')

        # Handle ms_levels default
        if ms_levels is None:
            ms_levels = self.getMSLevels() if num_spectra > 0 else []
        elif not ms_levels:
            ms_levels = self.getMSLevels() if num_spectra > 0 else []

        result = {}

        # Export spectra if requested
        if data in ('spectra', 'both'):
            spectra_table = self._build_spectra_arrow_table(
                format=format,
                columns=columns,
                ms_levels=ms_levels,
                min_rt=min_rt, max_rt=max_rt,
                min_mz=min_mz, max_mz=max_mz,
                include_precursor_info=include_precursor_info,
                include_ion_mobility=include_ion_mobility,
                pa=pa
            )
            result['spectra'] = spectra_table

        # Export chromatograms if requested
        if data in ('chromatograms', 'both'):
            chrom_table = self._build_chromatograms_arrow_table(
                format=format,
                columns=columns,
                min_rt=min_rt, max_rt=max_rt,
                pa=pa
            )
            result['chromatograms'] = chrom_table

        # Return single table or dict based on data parameter
        if data == 'both':
            return result
        elif data == 'spectra':
            return result['spectra']
        else:
            return result['chromatograms']

    def _build_spectra_arrow_table(self, format, columns, ms_levels, min_rt, max_rt,
                                    min_mz, max_mz, include_precursor_info,
                                    include_ion_mobility, pa):
        """Internal method to build Arrow table for spectra."""
        data_dict = {}

        if format == 'long':
            # Long format: one row per peak
            all_mz = []
            all_intensity = []
            all_rt = []
            all_spectrum_index = []
            all_ms_level = []
            all_native_id = []
            all_precursor_mz = []
            all_precursor_charge = []
            all_precursor_intensity = []
            all_isolation_lower = []
            all_isolation_upper = []
            all_ion_mobility = []

            for spec_idx, spec in enumerate(self):
                ms_level = spec.getMSLevel()
                if ms_level not in ms_levels:
                    continue

                rt = spec.getRT()
                if rt < min_rt or rt > max_rt:
                    continue

                mzs, intensities = spec.get_peaks()
                original_len = len(mzs)
                mz_mask = None

                # Apply m/z filter
                if len(mzs) > 0:
                    mz_mask = (mzs >= min_mz) & (mzs <= max_mz)
                    mzs = mzs[mz_mask]
                    intensities = intensities[mz_mask]

                n_peaks = len(mzs)
                if n_peaks == 0:
                    continue

                all_mz.append(mzs)
                all_intensity.append(intensities)
                all_rt.append(np.full(n_peaks, rt, dtype=np.float32))
                all_spectrum_index.append(np.full(n_peaks, spec_idx, dtype=np.uint32))
                all_ms_level.append(np.full(n_peaks, ms_level, dtype=np.uint8))
                all_native_id.append(np.full(n_peaks, spec.getNativeID(), dtype=object))

                if include_precursor_info:
                    precursors = spec.getPrecursors()
                    if precursors:
                        prec = precursors[0]
                        all_precursor_mz.extend([prec.getMZ()] * n_peaks)
                        all_precursor_charge.extend([prec.getCharge()] * n_peaks)
                        all_precursor_intensity.extend([prec.getIntensity()] * n_peaks)
                        all_isolation_lower.extend([prec.getIsolationWindowLowerOffset()] * n_peaks)
                        all_isolation_upper.extend([prec.getIsolationWindowUpperOffset()] * n_peaks)
                    else:
                        # Use None for proper Arrow null representation
                        all_precursor_mz.extend([None] * n_peaks)
                        all_precursor_charge.extend([None] * n_peaks)
                        all_precursor_intensity.extend([None] * n_peaks)
                        all_isolation_lower.extend([None] * n_peaks)
                        all_isolation_upper.extend([None] * n_peaks)

                if include_ion_mobility:
                    if spec.containsIMData():
                        im_idx, _ = spec.getIMData()
                        im_arr = spec.getFloatDataArrays()[im_idx]
                        im_data = [im_arr[i] for i in range(len(im_arr))]
                        # Apply same m/z filter mask if it was created
                        if mz_mask is not None and len(im_data) == original_len:
                            im_data = [im_data[i] for i in range(len(mz_mask)) if mz_mask[i]]
                        if len(im_data) == n_peaks:
                            all_ion_mobility.extend(im_data)
                        else:
                            all_ion_mobility.extend([None] * n_peaks)
                    else:
                        # Use None for proper Arrow null representation
                        all_ion_mobility.extend([None] * n_peaks)

            # Build data dict - order matches C++ ArrowExport implementation
            data_dict['mz'] = np.concatenate(all_mz) if all_mz else np.array([], dtype=np.float64)
            data_dict['intensity'] = np.concatenate(all_intensity) if all_intensity else np.array([], dtype=np.float32)
            data_dict['rt'] = np.concatenate(all_rt) if all_rt else np.array([], dtype=np.float32)

            # ion_mobility comes after rt (before spectrum_index) to match C++
            if include_ion_mobility:
                data_dict['ion_mobility'] = pa.array(all_ion_mobility, type=pa.float32())

            data_dict['spectrum_index'] = np.concatenate(all_spectrum_index) if all_spectrum_index else np.array([], dtype=np.uint32)
            data_dict['ms_level'] = np.concatenate(all_ms_level) if all_ms_level else np.array([], dtype=np.uint8)
            data_dict['native_id'] = np.concatenate(all_native_id) if all_native_id else np.array([], dtype=object)

            if include_precursor_info:
                # Use pa.array with explicit types for proper null handling
                data_dict['precursor_mz'] = pa.array(all_precursor_mz, type=pa.float64())
                data_dict['precursor_charge'] = pa.array(all_precursor_charge, type=pa.int16())
                data_dict['precursor_intensity'] = pa.array(all_precursor_intensity, type=pa.float32())
                data_dict['isolation_lower'] = pa.array(all_isolation_lower, type=pa.float64())
                data_dict['isolation_upper'] = pa.array(all_isolation_upper, type=pa.float64())

        else:
            # Semi-wide format: one row per spectrum, arrays for m/z and intensity
            all_mz = []
            all_intensity = []
            all_rt = []
            all_spectrum_index = []
            all_ms_level = []
            all_native_id = []
            all_precursor_mz = []
            all_precursor_charge = []
            all_precursor_intensity = []
            all_isolation_lower = []
            all_isolation_upper = []
            all_ion_mobility = []

            for spec_idx, spec in enumerate(self):
                ms_level = spec.getMSLevel()
                if ms_level not in ms_levels:
                    continue

                rt = spec.getRT()
                if rt < min_rt or rt > max_rt:
                    continue

                mzs, intensities = spec.get_peaks()
                original_len = len(mzs)
                mz_mask = None

                # Apply m/z filter
                if len(mzs) > 0:
                    mz_mask = (mzs >= min_mz) & (mzs <= max_mz)
                    mzs = mzs[mz_mask]
                    intensities = intensities[mz_mask]

                all_mz.append(mzs.tolist())
                all_intensity.append(intensities.tolist())
                all_rt.append(rt)
                all_spectrum_index.append(spec_idx)
                all_ms_level.append(ms_level)
                all_native_id.append(spec.getNativeID())

                if include_precursor_info:
                    precursors = spec.getPrecursors()
                    if precursors:
                        prec = precursors[0]
                        all_precursor_mz.append(prec.getMZ())
                        all_precursor_charge.append(prec.getCharge())
                        all_precursor_intensity.append(prec.getIntensity())
                        all_isolation_lower.append(prec.getIsolationWindowLowerOffset())
                        all_isolation_upper.append(prec.getIsolationWindowUpperOffset())
                    else:
                        # Use None for proper Arrow null representation
                        all_precursor_mz.append(None)
                        all_precursor_charge.append(None)
                        all_precursor_intensity.append(None)
                        all_isolation_lower.append(None)
                        all_isolation_upper.append(None)

                if include_ion_mobility:
                    if spec.containsIMData():
                        im_idx, _ = spec.getIMData()
                        im_arr = spec.getFloatDataArrays()[im_idx]
                        im_data = [im_arr[i] for i in range(len(im_arr))]
                        # Apply same m/z filter mask if it was created
                        if mz_mask is not None and len(im_data) == original_len:
                            im_data = [im_data[i] for i in range(len(mz_mask)) if mz_mask[i]]
                        all_ion_mobility.append(im_data)
                    else:
                        # Use None for null list (spectrum has no IM data)
                        all_ion_mobility.append(None)

            # Build data dict with list columns - order matches C++ ArrowExport
            data_dict['mz'] = all_mz
            data_dict['intensity'] = all_intensity
            data_dict['rt'] = np.array(all_rt, dtype=np.float32)

            # ion_mobility comes after rt (before spectrum_index) to match C++
            if include_ion_mobility:
                data_dict['ion_mobility'] = all_ion_mobility

            data_dict['spectrum_index'] = np.array(all_spectrum_index, dtype=np.uint32)
            data_dict['ms_level'] = np.array(all_ms_level, dtype=np.uint8)
            data_dict['native_id'] = all_native_id

            if include_precursor_info:
                # Use pa.array with explicit types for proper null handling
                data_dict['precursor_mz'] = pa.array(all_precursor_mz, type=pa.float64())
                data_dict['precursor_charge'] = pa.array(all_precursor_charge, type=pa.int16())
                data_dict['precursor_intensity'] = pa.array(all_precursor_intensity, type=pa.float32())
                data_dict['isolation_lower'] = pa.array(all_isolation_lower, type=pa.float64())
                data_dict['isolation_upper'] = pa.array(all_isolation_upper, type=pa.float64())

        # Filter columns if requested
        if columns is not None:
            data_dict = {k: v for k, v in data_dict.items() if k in columns}

        return pa.Table.from_pydict(data_dict)

    def _build_chromatograms_arrow_table(self, format, columns, min_rt, max_rt, pa):
        """Internal method to build Arrow table for chromatograms."""
        data_dict = {}

        if format == 'long':
            # Long format: one row per data point
            all_rt = []
            all_intensity = []
            all_chrom_index = []
            all_native_id = []
            all_precursor_mz = []
            all_product_mz = []

            for chrom_idx in range(self.getNrChromatograms()):
                chrom = self.getChromatogram(chrom_idx)
                rts, intensities = chrom.get_peaks()

                # Apply RT filter
                if len(rts) > 0:
                    mask = (rts >= min_rt) & (rts <= max_rt)
                    rts = rts[mask]
                    intensities = intensities[mask]

                n_points = len(rts)
                if n_points == 0:
                    continue

                all_rt.append(rts)
                all_intensity.append(intensities)
                all_chrom_index.append(np.full(n_points, chrom_idx, dtype=np.uint32))
                all_native_id.append(np.full(n_points, chrom.getNativeID(), dtype=object))
                all_precursor_mz.append(np.full(n_points, chrom.getPrecursor().getMZ(), dtype=np.float64))
                all_product_mz.append(np.full(n_points, chrom.getProduct().getMZ(), dtype=np.float64))

            # Build data dict
            data_dict['rt'] = np.concatenate(all_rt) if all_rt else np.array([], dtype=np.float64)
            data_dict['intensity'] = np.concatenate(all_intensity) if all_intensity else np.array([], dtype=np.float32)
            data_dict['chromatogram_index'] = np.concatenate(all_chrom_index) if all_chrom_index else np.array([], dtype=np.uint32)
            data_dict['native_id'] = np.concatenate(all_native_id) if all_native_id else np.array([], dtype=object)
            data_dict['precursor_mz'] = np.concatenate(all_precursor_mz) if all_precursor_mz else np.array([], dtype=np.float64)
            data_dict['product_mz'] = np.concatenate(all_product_mz) if all_product_mz else np.array([], dtype=np.float64)

        else:
            # Semi-wide format: one row per chromatogram
            all_rt = []
            all_intensity = []
            all_chrom_index = []
            all_native_id = []
            all_precursor_mz = []
            all_product_mz = []

            for chrom_idx in range(self.getNrChromatograms()):
                chrom = self.getChromatogram(chrom_idx)
                rts, intensities = chrom.get_peaks()

                # Apply RT filter
                if len(rts) > 0:
                    mask = (rts >= min_rt) & (rts <= max_rt)
                    rts = rts[mask]
                    intensities = intensities[mask]

                all_rt.append(rts.tolist())
                all_intensity.append(intensities.tolist())
                all_chrom_index.append(chrom_idx)
                all_native_id.append(chrom.getNativeID())
                all_precursor_mz.append(chrom.getPrecursor().getMZ())
                all_product_mz.append(chrom.getProduct().getMZ())

            # Build data dict with list columns
            data_dict['rt'] = all_rt
            data_dict['intensity'] = all_intensity
            data_dict['chromatogram_index'] = np.array(all_chrom_index, dtype=np.uint32)
            data_dict['native_id'] = all_native_id
            data_dict['precursor_mz'] = np.array(all_precursor_mz, dtype=np.float64)
            data_dict['product_mz'] = np.array(all_product_mz, dtype=np.float64)

        # Filter columns if requested
        if columns is not None:
            data_dict = {k: v for k, v in data_dict.items() if k in columns}

        return pa.Table.from_pydict(data_dict)

    def get_ion_df(self):
        """
        get_ion_df(self: MSExperiment) -> pd.DataFrame

        Generates a pandas DataFrame with all MS1 peaks and their ion mobility in the MSExperiment.

        Only MS level 1 spectra are exported.

        :return: Feature information stored in a DataFrame
        :rtype: pd.DataFrame

        :raises ImportError: If pandas is not installed
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for get_ion_df(). "
                "Please install it with: pip install pandas"
            )
        cols = ["rt", "mz", "intensity", "ion_mobility"]
        self.updateRanges()
        spectraarrs2d = self.get2DPeakDataIMLong(self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ(), 1)
        return pd.DataFrame(dict(zip(cols, spectraarrs2d)))

    def get_massql_df(self, ion_mobility=False):
        """
        get_massql_df(self: MSExperiment, ion_mobility: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]

        Exports data from MSExperiment to pandas DataFrames to be used with MassQL.

        The Python module massql allows queries in mass spectrometry data (MS1 and MS2
        data frames) in a SQL like fashion (https://github.com/mwang87/MassQueryLanguage).

        Both dataframes contain the columns:
        'i': intensity of a peak
        'i_norm': intensity normalized by the maximun intensity in the spectrum
        'i_tic_norm': intensity normalized by the sum of intensities (TIC) in the spectrum
        'mz': mass to charge of a peak
        'scan': number of the spectrum
        'rt': retention time of the spectrum
        'polarity': ion mode of the spectrum as integer value (positive: 1, negative: 2)
        'ion': the ionic mobility of a peak if ion parameter is True

        The MS2 dataframe contains additional columns:
        'precmz': mass to charge of the precursor ion
        'ms1scan': number of the corresponding MS1 spectrum
        'charge': charge of the precursor ion

        :param ion_mobility: If True, returns the ion mobility of the peaks.
        :type ion_mobility: bool

        :return: ms1_df - peak data of MS1 spectra, ms2_df - peak data of MS2 spectra with precursor information
        :rtype: Tuple[pd.DataFrame, pd.DataFrame]

        :raises ImportError: If pandas is not installed
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for get_massql_df(). "
                "Please install it with: pip install pandas"
            )
        from . import IonSource as _IonSource
        self.updateRanges()

        def _get_polarity(spec):
            '''Returns polarity as an integer value for the massql dataframe.

            According to massql positive polarity is represented by 1 and negative by 2.

            Parameters:
            spec (MSSpectrum): the spectrum to extract polarity

            Returns:
            int: polarity as int value according to massql specification
            '''
            polarity = spec.getInstrumentSettings().getPolarity()
            if polarity == _IonSource.Polarity.POLNULL:
                return 0
            elif polarity == _IonSource.Polarity.POSITIVE:
                return 1
            elif polarity == _IonSource.Polarity.NEGATIVE:
                return 2

        def _get_spec_arrays(mslevel):
            '''Get spectrum data as a matrix.

            Generator yields peak data from each spectrum (with specified MS level) as a numpy.ndarray.
            Normalized intensity values are calculated and the placeholder values replaced. For 'i_norm' and
            'i_tic_norm' the intensity values are divided by the maximum intensity value in the spectrum and
            the sum of intensity values, respectively.

            Parameters:
            mslevel (int): only spectra with the given MS level will be considered

            Yields:
            np.ndarray: 2D array with peak data (rows) from each spectrum
            '''
            for scan_num, spec in enumerate(self):
                if spec.getMSLevel() == mslevel:
                    mz, inty = spec.get_peaks()
                    # Safe division: handle empty spectra or all-zero intensities
                    max_inty = np.amax(inty, initial=0)
                    sum_inty = np.sum(inty)
                    i_norm = np.zeros_like(inty) if max_inty == 0 else inty / max_inty
                    i_tic_norm = np.zeros_like(inty) if sum_inty == 0 else inty / sum_inty
                    # data for both DataFrames: i, i_norm, i_tic_norm, mz, scan, rt, polarity
                    data = (inty, i_norm, i_tic_norm, mz, scan_num + 1, spec.getRT()/60, _get_polarity(spec))
                    cols = 7
                    if mslevel == 2:
                        cols = 10
                        # data for MS2 only: precmz, ms1scan, charge
                        # set fallback values if no precursor is annotated (-1)
                        if spec.getPrecursors():
                            data += (spec.getPrecursors()[0].getMZ(), self.getPrecursorSpectrum(scan_num)+1, spec.getPrecursors()[0].getCharge())
                        else:
                            data += (-1, -1, -1)
                    # create empty ndarr with shape according to MS level
                    ndarr = np.empty(shape=(spec.size(), cols))
                    # set column values
                    for i in range(cols):
                        ndarr[:,i] = data[i]
                    yield ndarr

        def _get_ion_spec_arrays(mslevel):
            '''Get spectrum data as a matrix.

            Generator yields peak data from each spectrum (with specified MS level) as a numpy.ndarray.
            Normalized intensity values are calculated and the placeholder values replaced. For 'i_norm' and
            'i_tic_norm' the intensity values are divided by the maximum intensity value in the spectrum and
            the sum of intensity values, respectively.

            Parameters:
            mslevel (int): only spectra with the given MS level will be considered

            Yields:
            np.ndarray: 2D array with peak data (rows) from each spectrum
            '''
            for scan_num, spec in enumerate(self):
                if spec.getMSLevel() == mslevel:
                    mz, inty = spec.get_peaks()
                    ion_array_idx, ion_unit = spec.getIMData()
                    ion_data_arr = spec.getFloatDataArrays()[ion_array_idx]
                    ion_data = ion_data_arr.get_data()

                    # Safe division: handle empty spectra or all-zero intensities
                    max_inty = np.amax(inty, initial=0)
                    sum_inty = np.sum(inty)
                    i_norm = np.zeros_like(inty) if max_inty == 0 else inty / max_inty
                    i_tic_norm = np.zeros_like(inty) if sum_inty == 0 else inty / sum_inty
                    # data for both DataFrames: i, i_norm, i_tic_norm, mz, scan, rt, polarity
                    data = (inty, i_norm, i_tic_norm, mz, scan_num + 1, spec.getRT()/60, _get_polarity(spec), ion_data)
                    cols = 8
                    if mslevel == 2:
                        cols = 11
                        # data for MS2 only: precmz, ms1scan, charge
                        # set fallback values if no precursor is annotated (-1)
                        if spec.getPrecursors():
                            data += (spec.getPrecursors()[0].getMZ(), self.getPrecursorSpectrum(scan_num)+1, spec.getPrecursors()[0].getCharge())
                        else:
                            data += (-1, -1, -1)
                    # create empty ndarr with shape according to MS level
                    ndarr = np.empty(shape=(spec.size(), cols))
                    # set column values
                    for i in range(cols):
                        ndarr[:,i] = data[i]
                    yield ndarr

        # create DataFrame for MS1 and MS2 with according column names and data types
        # if there are no spectra of given MS level return an empty DataFrame
        dtypes = {'i': 'float32', 'i_norm': 'float32', 'i_tic_norm': 'float32', 'mz': 'float64', 'scan': 'int32', 'rt': 'float32', 'polarity': 'int32'}
        if ion_mobility:
            dtypes = dict(dtypes, **{"ion_mobility": "float32"})

        if 1 in self.getMSLevels():
            spec_arrays = _get_spec_arrays(1) if not ion_mobility else _get_ion_spec_arrays(1)
            ms1_df = pd.DataFrame(np.concatenate(list(spec_arrays), axis=0), columns=dtypes.keys()).astype(dtypes)
        else:
            ms1_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

        dtypes = dict(dtypes, **{'precmz': 'float64', 'ms1scan': 'int32', 'charge': 'int32'})
        if 2 in self.getMSLevels():
            spec_arrays = _get_spec_arrays(2) if not ion_mobility else _get_ion_spec_arrays(2)
            ms2_df = pd.DataFrame(np.concatenate(list(spec_arrays), axis=0), columns=dtypes.keys()).astype(dtypes)
        else:
            ms2_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

        return ms1_df, ms2_df
