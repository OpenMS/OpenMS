

    def create(self):
        cdef IsotopeMarker res = IsotopeMarker.__new__(IsotopeMarker)
        res.inst = shared_ptr[_IsotopeMarker]( <_IsotopeMarker *> self.inst.get().create() )
        return res
