

    def transform(self, bytes path, transformer):
        cdef _String path_string = _String(<char *>path)
        assert hasattr(transformer, "consumeSpectrum"), "expected method consumeSpectrum"
        assert hasattr(transformer, "consumeChromatogram"), "expected method consumeChromatogram"
        assert hasattr(transformer, "setExpectedSize"), "expected method setExpectedSize"
        assert hasattr(transformer, "setExperimentalSettings"), "expected method setExperimentalSettings"
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum_mzml,
                                             _wrap_MSChromatogram_mzml,
                                             _wrap_ExperimentalSettings_mzml)

        try:
            self.inst.get().transform(path_string, consumer)
        finally:
            del consumer


cdef _wrap_MSSpectrum_mzml(const _MSSpectrum & _spec):
    cdef MSSpectrum spec = MSSpectrum.__new__(MSSpectrum)
    spec.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(_spec))
    return spec


cdef _wrap_MSChromatogram_mzml(const _MSChromatogram & _chromo):
    cdef MSChromatogram chromo = MSChromatogram.__new__(MSChromatogram)
    chromo.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(_chromo))
    return chromo


cdef _wrap_ExperimentalSettings_mzml(const _ExperimentalSettings & _exp):
    cdef ExperimentalSettings exp = ExperimentalSettings.__new__(ExperimentalSettings)
    exp.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings(_exp))
    return exp
