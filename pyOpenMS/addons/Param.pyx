


    def asDict(Param self):


        cdef list keys = list()

        cdef _ParamIterator param_it = self.inst.get().begin()
        while param_it != self.inst.get().end():
            keys.append(param_it.getName().c_str())
            inc(param_it)

        result = dict()

        for k in keys:
            value = self.getValue(k)
            dt = value.valueType()
            if dt == DataType.STRING_VALUE:
                value = value.toString()
            elif dt == DataType.INT_VALUE:
                value = value.toInt()
            elif dt == DataType.DOUBLE_VALUE:
                value = value.toDouble()
            elif dt == DataType.DOUBLE_LIST:
                value = value.toDoubleList()
            elif dt == DataType.STRING_LIST:
                value = value.toStringList()
            elif dt == DataType.INT_LIST:
                value = value.toIntList()
            result[k] = value

        return result

    def updateFrom(self, dd):
        assert isinstance(dd, dict)

        for key, v in dd.items():
            tags = self.getTags(key)
            desc = self.getDescription(key)
            self.setValue(key, DataValue(v), desc, tags)




