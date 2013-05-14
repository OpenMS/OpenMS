

    def getSpectrumById(self,  id ):
        assert isinstance(id, (int, long)), 'arg id wrong type'
    
        _r = self.inst.get().getSpectrumById((<int>id))
        cdef shared_ptr[_Spectrum] py_result = _r
        spec = Spectrum()
        spec.inst = _r
        return spec

    def getChromatogramById(self,  id ):
        assert isinstance(id, (int, long)), 'arg id wrong type'
    
        _r = self.inst.get().getChromatogramById((<int>id))
        cdef shared_ptr[_Chromatogram] py_result = _r
        chrom = Chromatogram()
        chrom.inst = _r
        return chrom
