

    def getAAFrequencies(self, dict mmap):
        cdef Map[_String, size_t] c_mmap 
        for k,v in mmap.iteritems():
            c_mmap[ _String(<char *>k) ] = v # add <size_t> ?
        self.inst.get().getAAFrequencies(c_mmap)

    def __iter__(self):
        cdef int n = self.inst.get().size()
        cdef int i = 0

        cdef Residue py_result 
        while i < n:

            py_result = Residue.__new__(Residue)
            py_result.inst = shared_ptr[_Residue](new _Residue( deref(self.inst.get())[i] ))

            yield py_result

            i += 1

