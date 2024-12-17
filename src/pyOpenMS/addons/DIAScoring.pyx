

    def dia_by_ion_score(self, list spectrum , AASequence sequence, charge, RangeMobility im_range, float bseries_score , float yseries_score):
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, OSSpectrum) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(sequence, AASequence), 'arg sequence wrong type'
        assert isinstance(im_range, RangeMobility), 'arg sequence wrong type'
        assert isinstance(charge, (int, long)), 'arg charge wrong type'
        assert isinstance(bseries_score, float), 'arg bseries_score wrong type'
        assert isinstance(yseries_score, float), 'arg yseries_score wrong type'

        cdef libcpp_vector[shared_ptr[_OSSpectrum]] v1
        cdef OSSpectrum spectrum_rec
        for spectrum_rec in spectrum:
            v1.push_back(spectrum_rec.inst)

        cdef double input_bseries_score = (<double>bseries_score)
        cdef double input_yseries_score = (<double>yseries_score)
        self.inst.get().dia_by_ion_score(v1, (deref(sequence.inst.get())), (<int>charge), (deref(im_range.inst.get())), input_bseries_score, input_yseries_score)
        yseries_score = input_yseries_score
        bseries_score = input_bseries_score
        return (bseries_score,yseries_score)

