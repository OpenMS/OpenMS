from DPosition cimport DPosition1 as _DPosition1
from DPosition cimport DPosition2 as _DPosition2

cdef class DPosition1:

    cdef shared_ptr[_DPosition1] inst

    def __dealloc__(self):
        self.inst.reset()

    def __init__(self):
        self.inst = shared_ptr[_DPosition1](new _DPosition1())

    def __getitem__(self, ix):
        if ix != 0:
            raise IndexError("invalid index %d" % ix)
        return deref(self.inst.get())[0]


cdef class DPosition2:

    cdef shared_ptr[_DPosition2] inst

    def __dealloc__(self):
        self.inst.reset()

    def __init__(self):
        self.inst = shared_ptr[_DPosition2](new _DPosition2())

    def __getitem__(self, ix):
        if ix != 0 and ix != 1:
            raise IndexError("invalid index %d" % ix)
        return deref(self.inst.get())[ix]
