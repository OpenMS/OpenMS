

    def convertToSpectrumPtr(self, MSSpectrum spectrum):
        assert isinstance(spectrum, MSSpectrum), 'arg spec wrong type'

        _r = self.inst.get().convertToSpectrumPtr((deref(spectrum.inst.get())))
        # cdef shared_ptr[_OSSpectrum] py_result = _r
        spec = OSSpectrum()
        spec.inst = _r
        return spec

    def convertToChromatogramPtr(self, MSChromatogram chrom):
        assert isinstance(chrom, MSChromatogram), 'arg spec wrong type'

        _r = self.inst.get().convertToChromatogramPtr((deref(chrom.inst.get())))
        # cdef shared_ptr[_OSChromatogram] py_result = _r
        ch = OSChromatogram()
        ch.inst = _r
        return ch

