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
    
    def descriptions(self):
        return [self.getDescription(k) for k in self.keys()]

    def update(self, *a):
        """
        use cases:

           p.update(dict d)
           p.update(Param p)
           p.update(Param p, int flag)
        """

        cdef Param p
        cdef int flag
        cdef bool r

        if len(a) == 1 and isinstance(a[0], dict):
            dd, = a
            for key, v in dd.items():
                self[key] = v
            return True
        elif len(a) == 1 and isinstance(a[0], Param):
            p, = a
            r = self.inst.get().update(<_Param>deref(p.inst.get()))
            return r
        elif len(a) == 2 and isinstance(a[0], Param) and isinstance(a[1], int):
            p, flag = a
            r = self.inst.get().update(<_Param>deref(p.inst.get()), <bool> flag)
            return r
        else:
            raise Exception("can not handle parameters of type %s" % (map(type, a)))

    def get(self, bytes key, default=None):
        if self.exists(key):
            return self.getValue(key)
        return default

    def __getitem__(self, bytes key):
        return self.getValue(key)

    def __setitem__(self, bytes key, value):
        tags = self.getTags(key)
        desc = self.getDescription(key)
        self.setValue(key, value, desc, tags)
        
    def __str__(self):
        return str(list(zip([k.decode() for k in self.keys()], self.values(), self.descriptions())))





