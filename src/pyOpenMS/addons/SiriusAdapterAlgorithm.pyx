# two empty lines are important !


    def preprocessingSirius(self,  featureinfo , MSExperiment spectra , list v_fp , KDTreeFeatureMaps fp_map_kd , FeatureMapping_FeatureToMs2Indices feature_mapping ):
        """Cython signature: void preprocessingSirius(const String & featureinfo, MSExperiment spectra, libcpp_vector[FeatureMap] & v_fp, KDTreeFeatureMaps & fp_map_kd, FeatureMapping_FeatureToMs2Indices & feature_mapping)"""
        assert (isinstance(featureinfo, str) or isinstance(featureinfo, unicode) or isinstance(featureinfo, bytes) or isinstance(featureinfo, String)), 'arg featureinfo wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert isinstance(v_fp, list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in v_fp), 'arg v_fp wrong type'
        assert isinstance(fp_map_kd, KDTreeFeatureMaps), 'arg fp_map_kd wrong type'
        assert isinstance(feature_mapping, FeatureMapping_FeatureToMs2Indices), 'arg feature_mapping wrong type'
    
    
        cdef libcpp_vector[_FeatureMap] * v2 = new libcpp_vector[_FeatureMap]()
        cdef FeatureMap item2
        for item2 in v_fp:
            v2.push_back(deref(item2.inst.get()))
    
    
        self.inst.get().preprocessingSirius(deref((convString(featureinfo)).get()), (deref(spectra.inst.get())), deref(v2), (deref(fp_map_kd.inst.get())), (deref(feature_mapping.inst.get())))
