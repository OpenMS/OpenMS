

    def getSpectrumAlignment(self, list result, spec1, spec2):
        # TODO: use autowrap for this directly
        assert isinstance(spec1, (MSSpectrum))
        assert isinstance(spec2, (MSSpectrum))
        cdef _MSSpectrum[_Peak1D] * _new_spec1 = new _MSSpectrum[_Peak1D]()
        cdef _MSSpectrum[_Peak1D] * _new_spec2 = new _MSSpectrum[_Peak1D]()
        cdef _Peak1D _peak

        if True:
            for _peak in deref(<_MSSpectrum[_Peak1D] *>(<MSSpectrum>spec1).inst.get()):
                _new_spec1.push_back(<_Peak1D>_peak)

        if True:
            for _peak in deref(<_MSSpectrum[_Peak1D] *>(<MSSpectrum>spec2).inst.get()):
                _new_spec2.push_back(<_Peak1D>_peak)

        cdef libcpp_vector[libcpp_pair[Size, Size]] _result

        self.inst.get().getSpectrumAlignment(_result, deref(_new_spec1), deref(_new_spec2))

        result[:] = []

        cdef libcpp_vector[libcpp_pair[Size, Size]].iterator _it =  _result.begin()
        while _it != _result.end():
            result.append((<int>(deref(_it).first), <int>(deref(_it).second)))
            inc(_it)

        del _new_spec1
        del _new_spec2



