from libcpp.map cimport map as libcpp_map


    def getAAFrequencies(self, dict mmap):
        cdef libcpp_map[_String, size_t] c_mmap
        self.inst.get().getAAFrequencies(c_mmap)
        for k,v in mmap.iteritems():
            v = c_mmap[ _String(<char *>k) ]

    def __iter__(self):
        cdef int n = self.inst.get().size()
        cdef int i = 0

        cdef Residue py_result 
        while i < n:

            py_result = Residue.__new__(Residue)
            py_result.inst = shared_ptr[_Residue](new _Residue( deref(self.inst.get())[i] ))

            yield py_result

            i += 1

