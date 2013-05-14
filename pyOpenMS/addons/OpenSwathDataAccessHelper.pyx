

    def convertToSpectrumPtr(self, MSSpectrum spectrum):
        assert isinstance(spectrum, MSSpectrum), 'arg spec wrong type'


        _r = self.inst.get().convertToSpectrumPtr((deref(spectrum.inst.get())))
        cdef shared_ptr[_Spectrum] py_result = _r
        spec = Spectrum()
        spec.inst = _r
        return spec

