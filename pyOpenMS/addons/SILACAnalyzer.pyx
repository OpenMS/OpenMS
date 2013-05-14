

    def initialize(self, bytes selected_labels_ ,  charge_min_ ,  charge_max_ ,  missed_cleavages_ ,  isotopes_per_peptide_min_ ,  isotopes_per_peptide_max_ , float rt_threshold_ , float rt_min_ , float intensity_cutoff_ , float intensity_correlation_ , float model_deviation_ ,  allow_missing_peaks_, label_identifiers):
        assert isinstance(selected_labels_, bytes), 'arg selected_labels_ wrong type'
        assert isinstance(charge_min_, (int, long)), 'arg charge_min_ wrong type'
        assert isinstance(charge_max_, (int, long)), 'arg charge_max_ wrong type'
        assert isinstance(missed_cleavages_, (int, long)), 'arg missed_cleavages_ wrong type'
        assert isinstance(isotopes_per_peptide_min_, (int, long)), 'arg isotopes_per_peptide_min_ wrong type'
        assert isinstance(isotopes_per_peptide_max_, (int, long)), 'arg isotopes_per_peptide_max_ wrong type'
        assert isinstance(rt_threshold_, float), 'arg rt_threshold_ wrong type'
        assert isinstance(rt_min_, float), 'arg rt_min_ wrong type'
        assert isinstance(intensity_cutoff_, float), 'arg intensity_cutoff_ wrong type'
        assert isinstance(intensity_correlation_, float), 'arg intensity_correlation_ wrong type'
        assert isinstance(model_deviation_, float), 'arg model_deviation_ wrong type'
        assert isinstance(allow_missing_peaks_, (int, long)), 'arg allow_missing_peaks_ wrong type'

        cdef libcpp_map[_String, DoubleReal] clabels
        for k,v in label_identifiers.iteritems():
          assert isinstance(k, bytes), 'arg key in label_identifiers wrong type'
          assert isinstance(v, float), 'arg value in label_identifiers wrong type'
          clabels[ _String(<char *>k) ] = (<float>v)

        self.inst.get().initialize((_String(<char *>selected_labels_)),
                                   (<unsigned int>charge_min_),
                                   (<unsigned int>charge_max_),
                                   (<int>missed_cleavages_),
                                   (<unsigned int>isotopes_per_peptide_min_),
                                   (<unsigned int>isotopes_per_peptide_max_),
                                   (<float>rt_threshold_),
                                   (<float>rt_min_),
                                   (<float>intensity_cutoff_),
                                   (<float>intensity_correlation_),
                                   (<float>model_deviation_),
                                   (<bool>allow_missing_peaks_),
                                   clabels ) 

