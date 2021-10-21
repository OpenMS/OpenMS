

    # define ArrayWrapper as holding in a vector
    cdef class ArrayWrapper:
        cdef libcpp_vector[float] vec
        cdef Py_ssize_t shape[1]
        cdef Py_ssize_t strides[1]

        # constructor and destructor are fairly unimportant now since
        # vec will be destroyed automatically.

        cdef set_data(self, libcpp_vector[float]& data):
           self.vec.swap(data)

        # now implement the buffer protocol for the class
        # which makes it generally useful to anything that expects an array
        def __getbuffer__(self, Py_buffer *buffer, int flags):
            # relevant documentation http://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
            cdef Py_ssize_t itemsize = sizeof(self.vec[0])

            self.shape[0] = self.vec.size()
            self.strides[0] = sizeof(int)
            buffer.buf = <char *>&(self.vec[0])
            buffer.format = 'f'
            buffer.internal = NULL
            buffer.itemsize = itemsize
            buffer.len = self.vec.size() * itemsize   # product(shape) * itemsize
            buffer.ndim = 1
            buffer.obj = self
            buffer.readonly = 0
            buffer.shape = self.shape
            buffer.strides = self.strides
            buffer.suboffsets = NULL


    def get2DPeakDataLong(self, min_rt, max_rt, min_mz, max_mz):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty] get2DPeakDataLong(float min_rt, float max_rt, float min_mz, float max_mz)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[float] mz
        cdef libcpp_vector[float] inty
        exp_.get2DPeakData(min_rt, max_rt, min_mz, max_mz, rt, mz, inty)
       
        cdef ArrayWrapper rt_wrap
        cdef ArrayWrapper mz_wrap
        cdef ArrayWrapper inty_wrap
        rt_wrap.set_data(rt)
        mz_wrap.set_data(mz)
        inty_wrap.set_data(inty)

        return (np.asarray(rt_wrap), np.asarray(mz_wrap),  np.asarray(inty_wrap))

    def getMSLevels(self):
        """Cython signature: list[int] getMSLevels()"""
        cdef libcpp_vector[unsigned int] _r = self.inst.get().getMSLevels()
        cdef libcpp_vector[unsigned int].iterator it__r = _r.begin()
        cdef list result = []
        while it__r != _r.end():
            result.append(deref(it__r))
            inc(it__r)
        return result

    def getChromatogram(self,  id_ ):
        """Cython signature: MSChromatogram getChromatogram(size_t id_)"""
        assert isinstance(id_, (int, long)), 'arg id_ wrong type'
        assert id_ < self.getNrChromatograms(), 'Requested chromatogram %s does not exist, there are only %s chromatograms' % (id_, self.getNrChromatograms() )
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogram((<size_t>id_)))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result

    def getSpectrum(self,  id_ ):
        """Cython signature: MSSpectrum getSpectrum(size_t id_)"""
        assert isinstance(id_, (int, long)), 'arg id_ wrong type'
        assert id_ < self.getNrSpectra(), 'Requested spectrum %s does not exist, there are only %s spectra' % (id_, self.getNrSpectra() )
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrum((<size_t>id_)))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
