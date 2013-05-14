# two empty lines are important !


    def asDict(Param self):


        cdef list keys = list()

        cdef _ParamIterator param_it = self.inst.get().begin()
        while param_it != self.inst.get().end():
            keys.append(param_it.getName().c_str())
            inc(param_it)

        result = dict()

        for k in keys:
            value = self.getValue(k)
            result[k] = value

        return result

    def updateFrom(self, dd):
        assert isinstance(dd, dict)

        for key, v in dd.items():
            tags = self.getTags(key)
            desc = self.getDescription(key)
            self.setValue(key, v, desc, tags)




