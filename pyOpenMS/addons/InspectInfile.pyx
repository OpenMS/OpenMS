


    def getModifications(self):
        _r = self.inst.get().getModifications()
        py_result = dict()
        cdef Map[_String, libcpp_vector[_String]].iterator outer_it = _r.begin()
        cdef libcpp_vector[_String].iterator inner_it
        cdef String item_0
        cdef str inner_key
        cdef list inner_values
        while outer_it != _r.end():
           inner_key = deref(outer_it).first.c_str()
           inner_values = []
           inner_it = deref(outer_it).second.begin()
           while inner_it != deref(outer_it).second.end():
               # item_0 = CVTerm.__new__(CVTerm)
               # item_0.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it)))
               inner_values.append(  deref(inner_it).c_str() )
               inc(inner_it)
           py_result[inner_key] = inner_values
           inc(outer_it)
        return py_result

