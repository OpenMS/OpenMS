

    def transform(self, bytes path, transformer):

        #
        # the referenced functions _wrap_MSSpectrum, _wrap_MSChromatogram and
        # _wrap_ExperimentalSettings are declared in the MzMLFile.pyx file!
        #

        cdef _String path_string = _String(<char *>path)
        assert hasattr(transformer, "consumeSpectrum")
        assert hasattr(transformer, "consumeChromatogram")
        assert hasattr(transformer, "setExpectedSize")
        assert hasattr(transformer, "setExperimentalSettings")
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum_mzxml,
                                             _wrap_MSChromatogram_mzxml,
                                             _wrap_ExperimentalSettings_mzxml)

        try:
            self.inst.get().transform(path_string, consumer)
        finally:
            del consumer

cdef _wrap_MSSpectrum_mzxml(const _MSSpectrum & _spec):
    cdef MSSpectrum spec = MSSpectrum.__new__(MSSpectrum)
    spec.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(_spec))
    return spec


cdef _wrap_MSChromatogram_mzxml(const _MSChromatogram & _chromo):
    cdef MSChromatogram chromo = MSChromatogram.__new__(MSChromatogram)
    chromo.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(_chromo))
    return chromo


cdef _wrap_ExperimentalSettings_mzxml(const _ExperimentalSettings & _exp):
    cdef ExperimentalSettings exp = ExperimentalSettings.__new__(ExperimentalSettings)
    exp.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings(_exp))
    return exp
