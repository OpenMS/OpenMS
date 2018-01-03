



    def optimizeCalibrationCurves(self, dict components_concentrations):
  
        # Input types:
        # components_concentrations   :  libcpp_map[String, libcpp_vector[AQS_featureConcentration]]

        # dict -> pair ( list, list ) )
        assert isinstance(components_concentrations, dict) and all(isinstance(elem, AQS_featureConcentration) 
            for li in components_concentrations.values()
                for elem in li), 'arg proteins wrong type'

        cdef libcpp_map[_String, libcpp_vector[ _AQS_featureConcentration] ] c_comp
        # declaration for the loop 
        cdef libcpp_vector[_AQS_featureConcentration] * c_vec_inner = new libcpp_vector[_AQS_featureConcentration]()
        cdef AQS_featureConcentration i_item1
        for k,v in components_concentrations.iteritems():
 
            c_vec_inner.clear()
 
            vec = v
            assert isinstance(vec, list) and all(isinstance(li, AQS_featureConcentration) for li in vec), 'arg proteins wrong type'
            for i_item1 in vec:
               c_vec_inner.push_back(deref(i_item1.inst.get()))
 
            assert isinstance(k, bytes), 'arg key in label_identifiers wrong type'
            c_comp[ _String(<char *>k) ] = deref(c_vec_inner)
 
        #
        ## Make the function call
        # 
        self.inst.get().optimizeCalibrationCurves(c_comp)

        #
        ## Get the data back from C++
        #
        replace = dict()
        cdef libcpp_map[_String, libcpp_vector[_AQS_featureConcentration] ].iterator it_map = c_comp.begin()
        cdef libcpp_vector[_AQS_featureConcentration] another_c_vec_inner 
        cdef libcpp_vector[_AQS_featureConcentration].iterator it_aqs
        cdef AQS_featureConcentration item_py
        cdef PeptideIdentification item_py_result_pep
        while it_map != c_comp.end():

            another_c_vec_inner = deref(it_map).second

            replace_inner = []
            it_aqs = another_c_vec_inner.begin()
            while it_aqs != another_c_vec_inner.end():
                item_py = AQS_featureConcentration.__new__(AQS_featureConcentration)   
                item_py.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration(deref(it_aqs)))
                replace_inner.append(item_py)
                inc(it_aqs)

            replace[ <libcpp_string>deref(it_map).first ]  = [ replace_inner] 
            inc(it_map)

        components_concentrations.clear()
        components_concentrations.update(replace)

        del c_vec_inner
 
