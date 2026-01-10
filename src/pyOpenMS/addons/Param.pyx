# two empty lines are important !


    @staticmethod
    def from_dict(d):
        """
        from_dict(d: dict) -> Param

        Create a new Param object from a Python dictionary.

        Args:
            d: Dictionary with string or bytes keys and parameter values

        Returns:
            A new Param object populated with the dictionary contents

        Example:
            >>> p = Param.from_dict({"algorithm:threshold": 0.5, "debug": True})
        """
        p = Param()
        p.update(d)
        return p

    def _to_bytes_key(self, key):
        """Convert a string or bytes key to bytes."""
        if isinstance(key, str):
            return key.encode('utf-8')
        return key

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

        Note: Consider using to_dict() for PEP 8 compliant naming.
        """
        return dict(self.items())

    def to_dict(self):
        """
        to_dict(self: Param) -> dict

        Convert the Param object to a Python dictionary.

        This is a PEP 8 compliant alias for asDict().

        Returns:
            Dictionary with string keys and parameter values

        Example:
            >>> p = Param()
            >>> p.setValue("threshold", 0.5, "A threshold value")
            >>> d = p.to_dict()
            >>> print(d)
            {'threshold': 0.5}
        """
        return self.asDict()

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

    def get(self, key, default=None):
        """
        get(self: Param, key: Union[str, bytes], default: Any = None) -> Any

        Get parameter value by key, with optional default.

        Args:
            key: Parameter key (string or bytes)
            default: Value to return if key doesn't exist (default: None)

        Returns:
            The parameter value, or default if not found

        Example:
            >>> value = param.get("threshold", 0.5)
        """
        key = self._to_bytes_key(key)
        if self.exists(key):
            return self.getValue(key)
        return default

    def __getitem__(self, key):
        """
        __getitem__(self: Param, key: Union[str, bytes]) -> Any

        Get parameter value by key using dictionary notation.

        Args:
            key: Parameter key (string or bytes)

        Returns:
            The parameter value

        Raises:
            Exception: If the key doesn't exist

        Example:
            >>> value = param["threshold"]
        """
        key = self._to_bytes_key(key)
        return self.getValue(key)

    def __setitem__(self, key, value):
        """
        __setitem__(self: Param, key: Union[str, bytes], value: Any) -> None

        Set parameter value by key using dictionary notation.

        Args:
            key: Parameter key (string or bytes)
            value: The value to set

        Example:
            >>> param["threshold"] = 0.5
        """
        key = self._to_bytes_key(key)
        tags = self.getTags(key)
        desc = self.getDescription(key)
        self.setValue(key, value, desc, tags)

    def __iter__(self):
        """
        __iter__(self: Param) -> Iterator[bytes]

        Iterate over parameter keys.

        Returns:
            Iterator over parameter keys

        Example:
            >>> for key in param:
            ...     print(key, param[key])
        """
        return iter(self.keys())

    def __len__(self):
        """
        __len__(self: Param) -> int

        Return the number of parameters.

        Returns:
            Number of parameters

        Example:
            >>> print(len(param))
        """
        return self.size()

    def __contains__(self, key):
        """
        __contains__(self: Param, key: Union[str, bytes]) -> bool

        Check if a parameter key exists.

        Args:
            key: Parameter key (string or bytes)

        Returns:
            True if key exists, False otherwise

        Example:
            >>> if "threshold" in param:
            ...     print("Found!")
        """
        key = self._to_bytes_key(key)
        return self.exists(key)

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
