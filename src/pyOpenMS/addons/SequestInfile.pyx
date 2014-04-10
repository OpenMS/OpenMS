

    def getModifications(self):
        cdef libcpp_map[_String, libcpp_vector[_String] ] c_res 
        c_res = self.inst.get().getModifications()
        pass 
        # 
        ## Get the data back from C++
        #
        replace = dict()
        # Make sure each C++ and each Cython objects you use are declared here
        # (only pure Python objects need no declaration).
        cdef libcpp_map[_String, libcpp_vector[_String] ].iterator it_outer = c_res.begin()
        cdef libcpp_vector[_String] c_inner_vec
        cdef libcpp_vector[_String].iterator it_inner
        cdef _String myString
        cdef String py_String
        while it_outer != c_res.end():
            
            inner_vec = []
            c_inner_vec = deref(it_outer).second
            it_inner = c_inner_vec.begin()
            while it_inner != c_inner_vec.end():

                py_String = String.__new__(String)   
                py_String.inst = shared_ptr[_String](new _String(deref(it_inner)))
                inner_vec.append(py_String)
                inc(it_inner)

            replace[ <libcpp_string>deref(it_outer).first ] = inner_vec
            inc(it_outer)

        return replace



