# two empty lines are important !


    def asDict(self):
        return dict(self.items())

    def keys(self):
        keys = list()
        cdef _ParamIterator param_it = self.inst.get().begin()
        while param_it != self.inst.get().end():
            key = param_it.getName().c_str()
            keys.append(key)
            inc(param_it)
        return keys

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def values(self):
        return [self[k] for k in self.keys()]

    def update(self, *a):

        cdef Param p
        cdef int flag

        if len(a) == 1 and isinstance(a[0], dict):
            dd, = a
            for key, v in dd.items():
                self[key] = v
        elif len(a) == 1 and isinstance(a[0], Param):
            p, = a
            self.inst.get().update(<_Param>deref(p.inst.get()))
        elif len(a) == 2 and isinstance(a[0], Param) and isinstance(a[1], int):
            p, flag = a
            self.inst.get().update(<_Param>deref(p.inst.get()), <bool> flag)
        else:
            raise Exception("can not handle parameters of type %s" % (map(type, a)))


    def __getitem__(self, str key):
        return self.getValue(key)

    def __setitem__(self, str key, value):
        tags = self.getTags(key)
        desc = self.getDescription(key)
        self.setValue(key, value, desc, tags)





