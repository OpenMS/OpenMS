

    def get_mz_array(MSSpectrum self):
        """Cython signature: numpy_vector get_mz_array()
        
        Will return a numpy array corresponding
        to the mz values in the MSSpectrum.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef np.ndarray[np.float64_t, ndim=1] mzs
        mzs = np.empty( (n,), dtype=np.float64)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            mzs[i] = deref(it).getMZ()
            inc(it)
            i += 1

        return mzs


    def get_intensity_array(MSSpectrum self):
        """Cython signature: numpy_vector get_intensity_array()
        
        Will return a numpy array corresponding
        to the intensity values in the MSSpectrum.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.empty( (n,), dtype=np.float32)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return intensities

    def get_peaks(self):
        """Cython signature: numpy_vector, numpy_vector get_peaks()
        
        Will return a tuple of two numpy arrays (m/z, intensity) corresponding
        to the peaks in the MSSpectrum. Provides fast access to peaks.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef np.ndarray[np.float64_t, ndim=1] mzs
        mzs = np.empty( (n,), dtype=np.float64)
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.empty( (n,), dtype=np.float32)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            mzs[i] = deref(it).getMZ()
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return mzs, intensities

    def set_peaks(self, peaks):
        """Cython signature: set_peaks((numpy_vector, numpy_vector))
        
        Takes a tuple or list of two arrays (m/z, intensity) and populates the
        MSSpectrum. The arrays can be numpy arrays (faster).
        """

        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"

        mzs, intensities = peaks
        assert len(mzs) == len(intensities), "Input vectors for set_peaks need to have the same length (mz and intensity vector)"

        # Select which function to use for set_peaks:
        # If we have numpy arrays, it helps to use optimized functions
        if isinstance(mzs, np.ndarray) and isinstance(intensities, np.ndarray) and \
          mzs.dtype == np.float64 and intensities.dtype == np.float32 and \
          mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_df(mzs, intensities)
        elif isinstance(mzs, np.ndarray) and isinstance(intensities, np.ndarray) and \
          mzs.dtype == np.float64 and intensities.dtype == np.float64 and \
          mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_dd(mzs, intensities)
        else:
            self._set_peaks_orig(mzs, intensities)



    def _set_peaks_fast_dd(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):

        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(data_mz)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef double intensity
        cdef int N
        N = len(data_mz)

        for i in range(N):
            mz = data_mz[i]
            intensity = data_i[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()


    def _set_peaks_fast_df(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[float, ndim=1, mode="c"] data_i not None):

        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(data_mz)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef float intensity
        cdef int N
        N = len(data_mz)

        for i in range(N):
            mz = data_mz[i]
            intensity = data_i[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()


    def _set_peaks_orig(self, mzs, intensities):


        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(mzs)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef float intensity
        cdef int N
        N = len(mzs)

        for i in range(N):
            mz = mzs[i]
            intensity = intensities[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()

    def intensityInRange(self, float mzmin, float mzmax):

        cdef double I

        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef int N = spec_.size()

        I = 0.0
        for i in range(N):
                if deref(spec_)[i].getMZ() >= mzmin:
                    break

        cdef _Peak1D * p
        for j in range(i, N):
                p = address(deref(spec_)[j])
                if p.getMZ() > mzmax:
                    break
                I += p.getIntensity()

        return I

    def getIMData(self):

        cdef libcpp_pair[Size, _DriftTimeUnit] r = self.inst.get().getIMData()

        pos = r.first
        unit = <int>r.second

        return (pos, unit)

    @property
    def mz(self):
        """
        Access m/z values with support for indexed assignment.

        Returns an MzAccessor object that behaves like a numpy array but
        supports indexed assignment operations. This enables direct updates
        to m/z values using boolean masks, integer arrays, or slices.

        Returns
        -------
        MzAccessor
            An accessor object that supports numpy-style operations and
            indexed assignment. Can be converted to a regular numpy array
            using np.array(spec.mz).

        Examples
        --------
        >>> spec = MSSpectrum()
        >>> spec.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 1500.0]))
        >>> spec.mz
        MzAccessor([100., 200., 300.])
        >>> spec.mz[0]
        100.0
        >>> # Indexed assignment (new feature)
        >>> spec.mz[0] = 150.0
        >>> spec.mz[1:3] = [250.0, 350.0]
        >>> spec.mz[spec.mz > 200] = 999.0

        See Also
        --------
        intensity : Access intensity values
        get_peaks : Get both m/z and intensity arrays
        """
        return MzAccessor(self)

    @mz.setter
    def mz(self, values):
        """
        Set m/z values while keeping intensities unchanged.

        Parameters
        ----------
        values : array-like
            Array of m/z values with length equal to number of peaks.

        Raises
        ------
        ValueError
            If length of values doesn't match number of peaks.
        """
        values_np = np.asarray(values, dtype=np.float64)
        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef size_t n = spec_.size()

        if len(values_np) != n:
            raise ValueError(f"Length of m/z values ({len(values_np)}) must match number of peaks ({n})")

        # Get current intensities and set new peaks
        _, intensities = self.get_peaks()
        self.set_peaks((values_np, intensities))

    @property
    def intensity(self):
        """
        Access intensity values with support for indexed assignment.

        Returns an IntensityAccessor object that behaves like a numpy array but
        supports indexed assignment operations. This enables direct updates
        to intensity values using boolean masks, integer arrays, or slices.

        Returns
        -------
        IntensityAccessor
            An accessor object that supports numpy-style operations and
            indexed assignment. Can be converted to a regular numpy array
            using np.array(spec.intensity).

        Examples
        --------
        >>> spec = MSSpectrum()
        >>> spec.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 1500.0]))
        >>> spec.intensity
        IntensityAccessor([1000., 2000., 1500.], dtype=float32)
        >>> spec.intensity.max()
        2000.0
        >>> # Indexed assignment (new feature)
        >>> spec.intensity[0] = 5000.0
        >>> spec.intensity[1:3] = [6000.0, 7000.0]
        >>> spec.intensity[spec.intensity > 1500] = 9999.0

        See Also
        --------
        mz : Access m/z values
        get_peaks : Get both m/z and intensity arrays
        """
        return IntensityAccessor(self)

    @intensity.setter
    def intensity(self, values):
        """
        Set intensity values while keeping m/z values unchanged.

        Parameters
        ----------
        values : array-like
            Array of intensity values with length equal to number of peaks.

        Raises
        ------
        ValueError
            If length of values doesn't match number of peaks.
        """
        values_np = np.asarray(values, dtype=np.float32)
        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef size_t n = spec_.size()

        if len(values_np) != n:
            raise ValueError(f"Length of intensity values ({len(values_np)}) must match number of peaks ({n})")

        # Get current m/z values and set new peaks
        mzs, _ = self.get_peaks()
        self.set_peaks((mzs, values_np))


# MzAccessor class for indexed assignment support
cdef class MzAccessor:
    """
    Proxy class for m/z value access with indexed assignment support.

    This class enables numpy-style indexed assignment to m/z values while
    maintaining backward compatibility with existing code. It intercepts write
    operations through __setitem__ and delegates read operations to numpy arrays.

    Examples
    --------
    >>> spec.mz[0] = 150.0  # Single index
    >>> spec.mz[1:3] = [250.0, 350.0]  # Slice
    >>> spec.mz[spec.mz > 200] = 999.0  # Boolean indexing
    >>> spec.mz[[0, 2]] = [150.0, 350.0]  # Fancy indexing
    >>> mz_array = np.array(spec.mz)  # Convert to regular numpy array
    """
    cdef MSSpectrum parent_spec

    def __init__(self, MSSpectrum parent):
        """Initialize MzAccessor with parent MSSpectrum."""
        self.parent_spec = parent

    def __getitem__(self, key):
        """
        Get m/z values for the given key.

        Supports all numpy indexing types: integers, slices, boolean masks,
        and integer arrays. Delegates to numpy array for read operations.
        """
        mz_values = self.parent_spec.get_mz_array()
        return mz_values[key]

    def __setitem__(self, key, values):
        """
        Set m/z values for the given key.

        Intercepts write operations to update the underlying MSSpectrum.
        Handles boolean masks, integer arrays, and slices.
        """
        # Get current m/z values
        mz_values = self.parent_spec.get_mz_array().copy()
        n = len(mz_values)

        # Convert key to indices
        if isinstance(key, (int, np.integer)):
            # Handle negative indices
            idx = key if key >= 0 else n + key
            if idx < 0 or idx >= n:
                raise IndexError(
                    f"Index {key} out of bounds. Valid range is [-{n}, {n - 1}]"
                )
            indices = np.array([idx], dtype=np.intp)
            values_array = np.array([values], dtype=np.float64)
        elif isinstance(key, slice):
            indices = np.arange(n)[key]
            values_array = np.asarray(values, dtype=np.float64)
        else:
            key_array = np.asarray(key)
            if key_array.dtype == np.bool_:
                if len(key_array) != n:
                    raise IndexError(
                        f"Boolean mask length ({len(key_array)}) must match "
                        f"number of peaks ({n})"
                    )
                indices = np.where(key_array)[0]
            else:
                # Handle negative indices in fancy indexing
                indices = key_array.astype(np.intp)
                indices = np.where(indices < 0, n + indices, indices)
            values_array = np.asarray(values, dtype=np.float64)

        # Validate indices (after normalizing negative indices)
        if len(indices) > 0:
            if np.any(indices < 0) or np.any(indices >= n):
                raise IndexError(
                    f"Index out of bounds. Valid range is [0, {n - 1}]"
                )

        # Handle broadcasting of single value
        if values_array.ndim == 0 or (values_array.ndim == 1 and len(values_array) == 1):
            if values_array.ndim == 0:
                values_array = np.full(len(indices), values_array.item(), dtype=np.float64)
            else:
                values_array = np.full(len(indices), values_array[0], dtype=np.float64)
        elif len(values_array) != len(indices):
            raise ValueError(
                f"Length of values ({len(values_array)}) must match "
                f"number of selected indices ({len(indices)})"
            )

        # Update m/z values at selected indices
        mz_values[indices] = values_array

        # Set back to spectrum using existing setter
        self.parent_spec.mz = mz_values

    def __array__(self, dtype=None, copy=None):
        """Enable numpy.array() conversion."""
        arr = self.parent_spec.get_mz_array()
        if dtype is not None:
            arr = arr.astype(dtype)
        return arr

    def __len__(self):
        """Return number of peaks."""
        return self.parent_spec.size()

    # Arithmetic operators for numpy compatibility
    def __add__(self, other):
        return self.parent_spec.get_mz_array() + other

    def __radd__(self, other):
        return other + self.parent_spec.get_mz_array()

    def __sub__(self, other):
        return self.parent_spec.get_mz_array() - other

    def __rsub__(self, other):
        return other - self.parent_spec.get_mz_array()

    def __mul__(self, other):
        return self.parent_spec.get_mz_array() * other

    def __rmul__(self, other):
        return other * self.parent_spec.get_mz_array()

    def __truediv__(self, other):
        return self.parent_spec.get_mz_array() / other

    def __rtruediv__(self, other):
        return other / self.parent_spec.get_mz_array()

    # Comparison operators
    def __eq__(self, other):
        return self.parent_spec.get_mz_array() == other

    def __ne__(self, other):
        return self.parent_spec.get_mz_array() != other

    def __lt__(self, other):
        return self.parent_spec.get_mz_array() < other

    def __le__(self, other):
        return self.parent_spec.get_mz_array() <= other

    def __gt__(self, other):
        return self.parent_spec.get_mz_array() > other

    def __ge__(self, other):
        return self.parent_spec.get_mz_array() >= other

    def __iter__(self):
        """Iterate over m/z values."""
        return iter(self.parent_spec.get_mz_array())

    def __repr__(self):
        """String representation of MzAccessor."""
        return f"MzAccessor({self.parent_spec.get_mz_array()})"

    def __str__(self):
        """String conversion of MzAccessor."""
        return str(self.parent_spec.get_mz_array())


# IntensityAccessor class for indexed assignment support
cdef class IntensityAccessor:
    """
    Proxy class for intensity value access with indexed assignment support.

    This class enables numpy-style indexed assignment to intensity values while
    maintaining backward compatibility with existing code. It intercepts write
    operations through __setitem__ and delegates read operations to numpy arrays.

    Examples
    --------
    >>> spec.intensity[0] = 5000.0  # Single index
    >>> spec.intensity[1:3] = [6000.0, 7000.0]  # Slice
    >>> spec.intensity[spec.intensity > 1500] = 9999.0  # Boolean indexing
    >>> spec.intensity[[0, 2]] = [5000.0, 7000.0]  # Fancy indexing
    >>> int_array = np.array(spec.intensity)  # Convert to regular numpy array
    """
    cdef MSSpectrum parent_spec

    def __init__(self, MSSpectrum parent):
        """Initialize IntensityAccessor with parent MSSpectrum."""
        self.parent_spec = parent

    def __getitem__(self, key):
        """
        Get intensity values for the given key.

        Supports all numpy indexing types: integers, slices, boolean masks,
        and integer arrays. Delegates to numpy array for read operations.
        """
        int_values = self.parent_spec.get_intensity_array()
        return int_values[key]

    def __setitem__(self, key, values):
        """
        Set intensity values for the given key.

        Intercepts write operations to update the underlying MSSpectrum.
        Handles boolean masks, integer arrays, and slices.
        """
        # Get current intensity values
        int_values = self.parent_spec.get_intensity_array().copy()
        n = len(int_values)

        # Convert key to indices
        if isinstance(key, (int, np.integer)):
            # Handle negative indices
            idx = key if key >= 0 else n + key
            if idx < 0 or idx >= n:
                raise IndexError(
                    f"Index {key} out of bounds. Valid range is [-{n}, {n - 1}]"
                )
            indices = np.array([idx], dtype=np.intp)
            values_array = np.array([values], dtype=np.float32)
        elif isinstance(key, slice):
            indices = np.arange(n)[key]
            values_array = np.asarray(values, dtype=np.float32)
        else:
            key_array = np.asarray(key)
            if key_array.dtype == np.bool_:
                if len(key_array) != n:
                    raise IndexError(
                        f"Boolean mask length ({len(key_array)}) must match "
                        f"number of peaks ({n})"
                    )
                indices = np.where(key_array)[0]
            else:
                # Handle negative indices in fancy indexing
                indices = key_array.astype(np.intp)
                indices = np.where(indices < 0, n + indices, indices)
            values_array = np.asarray(values, dtype=np.float32)

        # Validate indices (after normalizing negative indices)
        if len(indices) > 0:
            if np.any(indices < 0) or np.any(indices >= n):
                raise IndexError(
                    f"Index out of bounds. Valid range is [0, {n - 1}]"
                )

        # Handle broadcasting of single value
        if values_array.ndim == 0 or (values_array.ndim == 1 and len(values_array) == 1):
            if values_array.ndim == 0:
                values_array = np.full(len(indices), values_array.item(), dtype=np.float32)
            else:
                values_array = np.full(len(indices), values_array[0], dtype=np.float32)
        elif len(values_array) != len(indices):
            raise ValueError(
                f"Length of values ({len(values_array)}) must match "
                f"number of selected indices ({len(indices)})"
            )

        # Update intensity values at selected indices
        int_values[indices] = values_array

        # Set back to spectrum using existing setter
        self.parent_spec.intensity = int_values

    def __array__(self, dtype=None, copy=None):
        """Enable numpy.array() conversion."""
        arr = self.parent_spec.get_intensity_array()
        if dtype is not None:
            arr = arr.astype(dtype)
        return arr

    def __len__(self):
        """Return number of peaks."""
        return self.parent_spec.size()

    # Arithmetic operators for numpy compatibility
    def __add__(self, other):
        return self.parent_spec.get_intensity_array() + other

    def __radd__(self, other):
        return other + self.parent_spec.get_intensity_array()

    def __sub__(self, other):
        return self.parent_spec.get_intensity_array() - other

    def __rsub__(self, other):
        return other - self.parent_spec.get_intensity_array()

    def __mul__(self, other):
        return self.parent_spec.get_intensity_array() * other

    def __rmul__(self, other):
        return other * self.parent_spec.get_intensity_array()

    def __truediv__(self, other):
        return self.parent_spec.get_intensity_array() / other

    def __rtruediv__(self, other):
        return other / self.parent_spec.get_intensity_array()

    # Comparison operators
    def __eq__(self, other):
        return self.parent_spec.get_intensity_array() == other

    def __ne__(self, other):
        return self.parent_spec.get_intensity_array() != other

    def __lt__(self, other):
        return self.parent_spec.get_intensity_array() < other

    def __le__(self, other):
        return self.parent_spec.get_intensity_array() <= other

    def __gt__(self, other):
        return self.parent_spec.get_intensity_array() > other

    def __ge__(self, other):
        return self.parent_spec.get_intensity_array() >= other

    def __iter__(self):
        """Iterate over intensity values."""
        return iter(self.parent_spec.get_intensity_array())

    def __repr__(self):
        """String representation of IntensityAccessor."""
        return f"IntensityAccessor({self.parent_spec.get_intensity_array()})"

    def __str__(self):
        """String conversion of IntensityAccessor."""
        return str(self.parent_spec.get_intensity_array())

    # Expose common numpy array methods
    def max(self):
        """Return maximum intensity value."""
        return self.parent_spec.get_intensity_array().max()

    def min(self):
        """Return minimum intensity value."""
        return self.parent_spec.get_intensity_array().min()

    def sum(self):
        """Return sum of intensity values."""
        return self.parent_spec.get_intensity_array().sum()

    def mean(self):
        """Return mean of intensity values."""
        return self.parent_spec.get_intensity_array().mean()

