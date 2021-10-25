


    def get2DPeakDataLong(self, min_rt, max_rt, min_mz, max_mz):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty] get2DPeakDataLong(float min_rt, float max_rt, float min_mz, float max_mz)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[float] mz
        cdef libcpp_vector[float] inty
        exp_.get2DPeakData(min_rt, max_rt, min_mz, max_mz, rt, mz, inty)
       
        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat mz_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat inty_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)
        mz_wrap.set_data(mz)
        inty_wrap.set_data(inty)

        return (np.asarray(rt_wrap), np.asarray(mz_wrap), np.asarray(inty_wrap))

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
