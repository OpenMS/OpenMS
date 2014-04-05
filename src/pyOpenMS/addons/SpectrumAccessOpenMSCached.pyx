

    def getSpectrumById(self,  id ):
        assert isinstance(id, (int, long)), 'arg id wrong type'
    
        _r = self.inst.get().getSpectrumById((<int>id))
        cdef shared_ptr[_OSSpectrum] py_result = _r
        spec = OSSpectrum()
        spec.inst = _r
        return spec

    def getChromatogramById(self,  id ):
        assert isinstance(id, (int, long)), 'arg id wrong type'
    
        _r = self.inst.get().getChromatogramById((<int>id))
        cdef shared_ptr[_OSChromatogram] py_result = _r
        chrom = OSChromatogram()
        chrom.inst = _r
        return chrom
