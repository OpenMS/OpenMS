

    def transform(self, bytes path, transformer):
        cdef _String path_string = _String(<char *>path)
        assert hasattr(transformer, "consumeSpectrum"), "expected method consumeSpectrum"
        assert hasattr(transformer, "consumeChromatogram"), "expected method consumeChromatogram"
        assert hasattr(transformer, "setExpectedSize"), "expected method setExpectedSize"
        assert hasattr(transformer, "setExperimentalSettings"), "expected method setExperimentalSettings"
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum,
                                             _wrap_MSChromatogram,
                                             _wrap_ExperimentalSettings)

        try:
            self.inst.get().transform(path_string, consumer)
        finally:
            del consumer


cdef _wrap_MSSpectrum(const _MSSpectrum[_Peak1D] & _spec):
    cdef MSSpectrum spec = MSSpectrum.__new__(MSSpectrum)
    spec.inst = shared_ptr[_MSSpectrum[_Peak1D]](new _MSSpectrum[_Peak1D](_spec))
    return spec


cdef _wrap_MSChromatogram(const _MSChromatogram[_ChromatogramPeak] & _chromo):
    cdef MSChromatogram chromo = MSChromatogram.__new__(MSChromatogram)
    chromo.inst = shared_ptr[_MSChromatogram[_ChromatogramPeak]](new _MSChromatogram[_ChromatogramPeak](_chromo))
    return chromo


cdef _wrap_ExperimentalSettings(const _ExperimentalSettings & _exp):
    cdef ExperimentalSettings exp = ExperimentalSettings.__new__(ExperimentalSettings)
    exp.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings(_exp))
    return exp
