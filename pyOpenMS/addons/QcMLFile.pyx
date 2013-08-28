


    def map2csv(self, dict csv_table, String separator):
        cdef libcpp_map[_String, libcpp_map[_String, _String] ] c_dict_outer
        cdef libcpp_map[_String, _String] c_dict_inner
        for k,v in csv_table.iteritems():
            c_dict_inner.clear()
            for k_i,v_i in v.iteritems():
                c_dict_inner[ _String(<char *>k_i) ] = _String(<char *>v_i)
            c_dict_outer[ _String(<char *>k) ] = c_dict_inner

        
        cdef _String _r = self.inst.get().map2csv(c_dict_outer, deref(separator.inst.get()) )
        py_result = _cast_const_away(<char*>_r.c_str())
        return py_result

