




    ## def __iter__(self):
    ##     _r = self.inst.get().getSequence()
    ##     cdef libcpp_vector[const _Ribonucleotide *].iterator it__r = _r.begin()
    ##     cdef Ribonucleotide item_py_result
    ##     while it__r != _r.end():
    ##        item_py_result = Ribonucleotide.__new__(Ribonucleotide)
    ##        item_py_result.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(deref(it__r))))
    ##        yield py_result
    ##        inc(it__r)

    def __iter__(self):

        for r in self.getSequence():
          yield r

