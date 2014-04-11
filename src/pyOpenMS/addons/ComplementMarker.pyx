

    def create(self):
        cdef ComplementMarker res = ComplementMarker.__new__(ComplementMarker)
        res.inst = shared_ptr[_ComplementMarker]( <_ComplementMarker *> self.inst.get().create() )
        return res
