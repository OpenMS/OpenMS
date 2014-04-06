

    def transform(self, bytes path, transformer):

        #
        # the refernced functions _wrap_MSSpectrum, _wrap_MSChromatogram and
        # _wrap_ExperimentalSettings are declared in MzMLFile.pyx next to this file !!!
        #

        cdef _String path_string = _String(<char *>path)
        assert hasattr(transformer, "consumeSpectrum")
        assert hasattr(transformer, "consumeChromatogram")
        assert hasattr(transformer, "setExpectedSize")
        assert hasattr(transformer, "setExperimentalSettings")
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum,
                                             _wrap_MSChromatogram,
                                             _wrap_ExperimentalSettings)

        try:
            self.inst.get().transform(path_string, consumer)
        finally:
            del consumer
