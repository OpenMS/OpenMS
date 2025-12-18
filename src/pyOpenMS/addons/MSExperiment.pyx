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

        Args:
            long_format (bool): If True, returns columns for long format.
                               If False, returns columns for compact format.

        Returns:
            list: List of column name strings.

        Example:
            >>> exp.get_df_columns(long_format=True)
            ['rt', 'mz', 'intensity', 'ms_level']

            >>> exp.get_df_columns(long_format=False)
            ['rt', 'ms_level', 'mz_array', 'intensity_array']
        """
        if long_format:
            return ['rt', 'mz', 'intensity', 'ms_level']
        else:
            return ['rt', 'ms_level', 'mz_array', 'intensity_array']

    def get_df(self, columns=None, ms_levels=None, long_format=False):
        """
        get_df(self: MSExperiment, columns: Optional[List[str]] = None, ms_levels: List[int] = [], long_format: bool = False) -> pd.DataFrame

        Generates a pandas DataFrame with all peaks in the MSExperiment

        Parameters:
        columns (list or None): List of column names to include. If None,
                               includes all columns. Use get_df_columns() to discover available columns.
        ms_levels (List[int]): Get only spectra with the given MS levels. Default is an empty list, which means all MS levels will be included.
        long_format (bool): set to True if you want to have a long/expanded/melted dataframe with one row per peak. Faster but
            replicated RT information. If False, returns rows in the style: rt, np.array(mz), np.array(int)

        Returns:
        pd.DataFrame: feature information stored in a DataFrame

        Example:
            >>> # Get all columns
            >>> df = exp.get_df()

            >>> # Discover available columns
            >>> print(exp.get_df_columns())

            >>> # Get only specific columns
            >>> df = exp.get_df(columns=['rt', 'mz', 'intensity'], long_format=True)
        """
        import pandas as pd
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

    def to_arrow(self, columns=None, ms_levels=None, long_format=False):
        """
        to_arrow(self: MSExperiment, columns: Optional[List[str]] = None, ms_levels: List[int] = [], long_format: bool = False) -> pa.Table

        Returns an Apache Arrow Table with all peaks in the MSExperiment.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all columns. Use get_df_columns() to discover available columns.
            ms_levels (List[int]): Get only spectra with the given MS levels. Default is an empty list,
                                  which means all MS levels will be included.
            long_format (bool): Set to True for a long/expanded table with one row per peak.
                              If False, returns rows with array columns (rt, ms_level, mz_array, intensity_array).

        Returns:
            pyarrow.Table: Arrow Table with peak data.

        Example:
            >>> # Get all columns in long format
            >>> table = exp.to_arrow(long_format=True)

            >>> # Convert to pandas (zero-copy with pandas 2.0+)
            >>> df = table.to_pandas()

            >>> # Convert to polars
            >>> import polars as pl
            >>> df = pl.from_arrow(table)
        """
        import pyarrow as pa
        if ms_levels is None:
            ms_levels = []
        self.updateRanges()
        if not ms_levels:
            ms_levels = self.getMSLevels()

        if long_format:
            # Build data dict for long format
            all_rts = []
            all_mzs = []
            all_intensities = []
            all_ms_levels = []

            for ms_level in ms_levels:
                rt_arr, mz_arr, inty_arr = self.get2DPeakDataLong(
                    self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ(), ms_level
                )
                all_rts.append(rt_arr)
                all_mzs.append(mz_arr)
                all_intensities.append(inty_arr)
                all_ms_levels.append(np.full(len(rt_arr), ms_level, dtype=np.int32))

            data_dict = {
                'rt': np.concatenate(all_rts) if all_rts else np.array([], dtype=np.float32),
                'mz': np.concatenate(all_mzs) if all_mzs else np.array([], dtype=np.float32),
                'intensity': np.concatenate(all_intensities) if all_intensities else np.array([], dtype=np.float32),
                'ms_level': np.concatenate(all_ms_levels) if all_ms_levels else np.array([], dtype=np.int32),
            }
        else:
            # Compact format with arrays
            rts = []
            ms_level_list = []
            mz_arrays = []
            intensity_arrays = []

            for spec in self:
                if spec.getMSLevel() in ms_levels:
                    rts.append(spec.getRT())
                    ms_level_list.append(spec.getMSLevel())
                    mz, inty = spec.get_peaks()
                    mz_arrays.append(mz)
                    intensity_arrays.append(inty)

            data_dict = {
                'rt': np.array(rts, dtype=np.float64),
                'ms_level': np.array(ms_level_list, dtype=np.int32),
                'mz_array': mz_arrays,  # List of arrays for Arrow list type
                'intensity_array': intensity_arrays,
            }

        # Filter columns if requested
        if columns is not None:
            data_dict = {k: v for k, v in data_dict.items() if k in columns}

        return pa.Table.from_pydict(data_dict)

    def get_ion_df(self):
        """
        get_ion_df(self: MSExperiment) -> pd.DataFrame

        Generates a pandas DataFrame with all peaks and the ion mobility in the MSExperiment

        Returns:
        pd.DataFrame: feature information stored in a DataFrame
        """
        import pandas as pd
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

        Parameters:
            ion_mobility (bool): if True, returns the ion mobility of the peaks.

        Returns:
        ms1_df (pd.DataFrame): peak data of MS1 spectra
        ms2_df (pd.DataFrame): peak data of MS2 spectra with precursor information
        """
        import pandas as pd
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
