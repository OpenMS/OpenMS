


    def transform(self, *args):
        if (len(args)==2):
             self._transform_1(*args)
        elif (len(args)==3):
             self._transform_2(*args)
        elif (len(args)==4):
             self._transform_3(*args)
        elif (len(args)==5):
             self._transform_4(*args)
        else:
           raise Exception('can not handle type of %s' % (args,))

        # void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *) nogil except + # wrap-ignore
        # void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, bool skip_full_count, bool skip_first_pass) nogil except +
        # void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, MSExperiment& e) nogil except + # wrap-ignore
        # void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, MSExperiment& e, bool skip_full_count, bool skip_first_pass) nogil except + # wrap-ignore

    def _transform_4(self, path, transformer, MSExperiment exp, bool skip_full_count, bool skip_first_pass):
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

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
            self.inst.get().transform(deref((convString(path)).get()), consumer, deref(exp.inst.get()), skip_full_count, skip_first_pass)
        finally:
            del consumer

    def _transform_2(self, path, transformer, MSExperiment exp):
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

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
            self.inst.get().transform(deref((convString(path)).get()), consumer, deref(exp.inst.get()) )
        finally:
            del consumer

    def _transform_3(self, path, transformer, bool skip_full_count, bool skip_first_pass):
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

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
            self.inst.get().transform(deref((convString(path)).get()), consumer, skip_full_count, skip_first_pass)
        finally:
            del consumer


    def _transform_1(self, path, transformer):
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

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
            self.inst.get().transform(deref((convString(path)).get()), consumer)
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
