# two empty lines are important !


    def initPluginParam(self, name, version):
        """
        initPluginParam(self: Param, name: Union[str, bytes], version: Union[str, bytes]) -> None
        
        Initialize parameters for a plugin with standard structure.
        """
        self.addSection(name, "Main section for Plugin '" + name + "'")
        self.setValue(name + ":version", version, "Version of Plugin '" + name + "'")
        self.addSection(name + ":1", "Instance '1' section for Plugin '" + name + "'")

        self.setValue(name + ":1:in", "", "The input file", [b"required", b"input file"])
        self.setValue(name + ":1:out", "", "The output file", [b"required", b"output file"])


    def asDict(self):
        """
        asDict(self: Param) -> dict
        
        Convert the Param object to a Python dictionary.
        """
        return dict(self.items())

    def keys(self):
        """
        keys(self: Param) -> List[bytes]
        
        Return a list of all parameter keys.
        """
        keys = list()
        cdef _ParamIterator param_it = self.inst.get().begin()
        while param_it != self.inst.get().end():
            key = param_it.getName().c_str()
            keys.append(key)
            inc(param_it)
        return keys

    def items(self):
        """
        items(self: Param) -> List[Tuple[bytes, Union[int, float, bytes, str, List[int], List[float], List[bytes]]]]
        
        Return a list of (key, value) tuples.
        """
        return [(k, self[k]) for k in self.keys()]

    def values(self):
        """
        values(self: Param) -> List[Union[int, float, bytes, str, List[int], List[float], List[bytes]]]
        
        Return a list of all parameter values.
        """
        return [self[k] for k in self.keys()]
    
    def descriptions(self):
        """
        descriptions(self: Param) -> List[bytes]
        
        Return a list of descriptions for all parameters.
        """
        return [self.getDescription(k) for k in self.keys()]

    def update(self, *a):
        """
        update(self: Param, *args: Union[dict, Param, Tuple[Param, int]]) -> bool
        
        Update parameter values from a dictionary or another Param object.
        
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
        """
        get(self: Param, key: bytes, default: Optional[Union[int, float, bytes, str, List[int], List[float], List[bytes]]] = None) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        
        Get parameter value by key, with optional default.
        """
        if self.exists(key):
            return self.getValue(key)
        return default

    def __getitem__(self, bytes key):
        """
        __getitem__(self: Param, key: bytes) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        
        Get parameter value by key using dictionary notation.
        """
        return self.getValue(key)

    def __setitem__(self, bytes key, value):
        """
        __setitem__(self: Param, key: bytes, value: Union[int, float, bytes, str, List[int], List[float], List[bytes]]) -> None
        
        Set parameter value by key using dictionary notation.
        """
        tags = self.getTags(key)
        desc = self.getDescription(key)
        self.setValue(key, value, desc, tags)
        
    def __str__(self):
        """
        __str__(self: Param) -> str
        
        Return a string representation of the Param object.
        """
        return str(list(zip([k.decode() for k in self.keys()], self.values(), self.descriptions())))

    def __repr__(self):
        """
        __repr__(self: Param) -> str
        
        Return a string representation of the Param object.
        """
        return self.__str__()
