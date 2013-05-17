

    def getAAFrequencies(self, dict mmap):
        cdef Map[_String, size_t] c_mmap 
        for k,v in mmap.iteritems():
            c_mmap[ _String(<char *>k) ] = v # add <size_t> ?
        self.inst.get().getAAFrequencies(c_mmap)
