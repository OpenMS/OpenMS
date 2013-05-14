# the following empty line is important

    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def getFileDescriptions(self):
        cdef FileDescriptions _r = self.inst.get().getFileDescriptions()
        py_result = dict()
        cdef FileDescriptions_iterator it__r = _r.begin()
        cdef FileDescription item_py_result
        while it__r != _r.end():
           item_py_result = FileDescription.__new__(FileDescription)
           item_py_result.inst = shared_ptr[_FileDescription](new _FileDescription((deref(it__r)).second))
           py_result[<unsigned long int>(deref(it__r).first)] = item_py_result
           inc(it__r)
        return py_result
    
    def setFileDescriptions(self, dict in_0 ):
        assert isinstance(in_0, dict) and all(isinstance(k, (int, long)) for k in in_0.keys()) and all(isinstance(v, FileDescription) for v in in_0.values()), 'arg in_0 wrong type'
        cdef FileDescriptions * v0 = new FileDescriptions()
        for key, value in in_0.items():
           deref(v0)[<unsigned long int> key] = deref((<FileDescription>value).inst.get())
        self.inst.get().setFileDescriptions(deref(v0))
        cdef replace_in_0 = dict()
        cdef FileDescriptions_iterator it_in_0 = v0.begin()
        cdef FileDescription item_in_0
        while it_in_0 != v0.end():
           item_in_0 = FileDescription.__new__(FileDescription)
           item_in_0.inst = shared_ptr[_FileDescription](new _FileDescription((deref(it_in_0)).second))
           replace_in_0[<unsigned long int> deref(it_in_0).first] = item_in_0
           inc(it_in_0)
        in_0.clear()
        in_0.update(replace_in_0)
        del v0
    
