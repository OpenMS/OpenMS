



    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_RibonucleotideDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_RibonucleotideDB](_getInstance_RibonucleotideDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass


    def getRibonucleotideAlternatives(self, bytes code ):
        """Cython signature: libcpp_pair[const Ribonucleotide *,const Ribonucleotide *] getRibonucleotideAlternatives(const libcpp_string & code)"""
        assert isinstance(code, bytes), 'arg code wrong type'
    
        _r = self.inst.get().getRibonucleotideAlternatives((<libcpp_string>code))
        cdef const _Ribonucleotide * out_ptr1 = _r.first
        cdef const _Ribonucleotide * out_ptr2 = _r.second

        cdef Ribonucleotide out1 = Ribonucleotide.__new__(Ribonucleotide)
        cdef Ribonucleotide out2 = Ribonucleotide.__new__(Ribonucleotide)

        out1.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(_r.first)))
        out2.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(_r.second)))

        cdef list py_result = [out1, out2]
        return py_result 
