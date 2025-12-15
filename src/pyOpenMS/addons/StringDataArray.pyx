



    def __getitem__(self,  in_0 ):
        """
        __getitem__(self: StringDataArray, in_0: int) -> bytes
        
        Get element at index in_0.
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
        assert in_0 >= 0, 'arg in_0 cannot be negative'

        cdef unsigned int _idx = (<int>in_0)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef _String _r = deref(self.inst.get())[(<int>in_0)]
        py_result = _cast_const_away(<char*>_r.c_str())
        return py_result

    def __setitem__(self, key, value):
        """
        __setitem__(self: StringDataArray, key: int, value: Union[str, bytes, unicode, String]) -> None
        
        Set element at index key to value.
        """
        assert isinstance(key, int), 'arg key wrong type'
        assert (isinstance(value, str) or isinstance(value, unicode) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
        assert key >= 0, 'arg key cannot be negative'

        cdef unsigned int _idx = (<int>key)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef shared_ptr[_String] _s = convString(value)
        deref(self.inst.get())[(<int>key)] = deref(_s.get())

    def __len__(self):
        """
        __len__(self: StringDataArray) -> int
        
        Return the number of elements in the array.
        """
        return self.inst.get().size()

    def __str__(self):
        """
        __str__(self: StringDataArray) -> str
        
        Return a string representation of the StringDataArray object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: StringDataArray) -> str
        
        Return a string representation of the StringDataArray object.

        Returns key properties in a readable format:
        StringDataArray(name='annotation', size=100)
        """
        cdef unsigned int n = self.inst.get().size()
        
        parts = []
        
        # Get name from MetaInfoDescription
        name = self.getName()
        if name:
            parts.append(f"name='{name}'")
        
        parts.append(f"size={n}")
        
        return f"StringDataArray({', '.join(parts)})"

    def get_data(self):
        """
        get_data(self: StringDataArray) -> List[bytes]
        
        Gets the data as a list of Python strings.

        This method creates a copy of the underlying data, so it's safe to use
        even after the original StringDataArray object is deleted or modified.

        Returns:
            list: A list of Python strings containing a copy of the data.

        Example usage:

        .. code-block:: python

            sda = pyopenms.StringDataArray()
            sda.push_back('a')
            sda.push_back('b')
            data = sda.get_data()  # ['a', 'b']
        """
        cdef _StringDataArray * sda_ = self.inst.get()
        cdef unsigned int n = sda_.size()

        if n == 0:
            return []

        result = []
        cdef unsigned int i
        cdef _String s
        for i in range(n):
            s = deref(sda_)[i]
            result.append(_cast_const_away(<char*>s.c_str()))

        return result

