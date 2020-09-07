

    def __call__(self, *args):
        assert len(args) in (1, 2), len(args)
        if len(args) == 1:
            return self._score_1(args[0])
        return self._score_2(args[0], args[1])

    def _score_1(self, spec1):
        assert isinstance(spec1, (MSSpectrum,))
        cdef _MSSpectrum * _new_spec1 = new _MSSpectrum()
        cdef _Peak1D _peak

        if True:
            for _peak in deref(<_MSSpectrum *>(<MSSpectrum>spec1).inst.get()):
                _new_spec1.push_back(<_Peak1D>_peak)

        cdef _SpectrumAlignmentScore * scorer = self.inst.get()
        cdef double score = <double>deref(scorer)(deref(_new_spec1))
        del _new_spec1
        return score

    def _score_2(self, spec1, spec2):
        assert isinstance(spec1, (MSSpectrum))
        assert isinstance(spec2, (MSSpectrum))
        cdef _MSSpectrum * _new_spec1 = new _MSSpectrum()
        cdef _MSSpectrum * _new_spec2 = new _MSSpectrum()
        cdef _Peak1D _peak

        if True:
            for _peak in deref(<_MSSpectrum *>(<MSSpectrum>spec1).inst.get()):
                _new_spec1.push_back(<_Peak1D>_peak)

        if True:
            for _peak in deref(<_MSSpectrum *>(<MSSpectrum>spec2).inst.get()):
                _new_spec2.push_back(<_Peak1D>_peak)

        cdef _SpectrumAlignmentScore * scorer = self.inst.get()
        cdef double score = <double>deref(scorer)(deref(_new_spec1), deref(_new_spec2))
        del _new_spec1
        del _new_spec2
        return score
