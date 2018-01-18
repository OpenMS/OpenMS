

    def dia_by_ion_score(self, OSSpectrum spectrum , AASequence sequence ,  charge , float bseries_score , float yseries_score ):
        assert isinstance(spectrum, OSSpectrum), 'arg spectrum wrong type'
        assert isinstance(sequence, AASequence), 'arg sequence wrong type'
        assert isinstance(charge, (int, long)), 'arg charge wrong type'
        assert isinstance(bseries_score, float), 'arg bseries_score wrong type'
        assert isinstance(yseries_score, float), 'arg yseries_score wrong type'
        cdef shared_ptr[_OSSpectrum] input_spectrum = spectrum.inst

        cdef double input_bseries_score = (<double>bseries_score)
        cdef double input_yseries_score = (<double>yseries_score)
        self.inst.get().dia_by_ion_score(input_spectrum, (deref(sequence.inst.get())), (<int>charge), input_bseries_score, input_yseries_score)
        yseries_score = input_yseries_score
        bseries_score = input_bseries_score
        return (bseries_score,yseries_score)


