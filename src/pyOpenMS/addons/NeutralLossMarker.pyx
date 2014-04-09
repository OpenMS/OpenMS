

    def create(self):
        cdef NeutralLossMarker res = NeutralLossMarker.__new__(NeutralLossMarker)
        res.inst = shared_ptr[_NeutralLossMarker]( <_NeutralLossMarker *> self.inst.get().create() )
        return res
