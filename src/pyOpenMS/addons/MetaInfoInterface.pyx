



    def getMetaValues(self):
        """Cython signature: dict getMetaValues()

        Returns all meta values in a Python dictionary.
        """
        # cdef _DataValue _c_value
        # cdef DataValue _value 
        # cdef int _type

        mmap = {}
        cdef libcpp_vector[_String] keys 
        cdef object py_result

        # get keys for iteration
        self.inst.get().getKeys(keys)

        cdef libcpp_vector[_String].iterator k_it = keys.begin()
        while k_it != keys.end():
            # easy approach: call Python fxn
            py_str = _cast_const_away(<char*>(deref(k_it)).c_str())
            py_result = self.getMetaValue(py_str)
            # # hard approach: do it ourselves
            # _c_value = self.inst.get().getMetaValue(deref(k_it))
            # py_str = _cast_const_away(<char*>(deref(k_it)).c_str())
            # _type = _c_value.valueType()
            # _value = DataValue.__new__(DataValue)
            # _value.inst = shared_ptr[_DataValue](new _DataValue(_c_value))
            # if _type == DataType.STRING_VALUE:
            #     py_result = _value.toString()
            # elif _type == DataType.INT_VALUE:
            #     py_result = _value.toInt()
            # elif _type == DataType.DOUBLE_VALUE:
            #     py_result = _value.toDouble()
            # elif _type == DataType.INT_LIST:
            #     py_result = _value.toIntList()
            # elif _type == DataType.DOUBLE_LIST:
            #     py_result = _value.toDoubleList()
            # elif _type == DataType.STRING_LIST:
            #     py_result = _value.toStringList()
            # elif _type == DataType.EMPTY_VALUE:
            #     py_result = None
            # else:
            #     raise Exception("DataValue instance has invalid value type %d" % _type)

            mmap[ py_str ] = py_result
            inc(k_it)

        return mmap

    def setMetaValues(self, dict mmap):
        """Cython signature: setMetaValues(dict values)

        Sets the meta values given in the Python dictionary.
        """
        self.inst.get().clearMetaInfo() # ensure its empty first
        for k, v in mmap.iteritems():
            self.inst.get().setMetaValue(deref((convString(k)).get()), deref(DataValue(v).inst.get()))
            # self.setMetaValue(k, v)






