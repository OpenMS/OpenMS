# This goes into the PXD file
from DPosition cimport DPosition1 as _DPosition1
from DPosition cimport DPosition2 as _DPosition2
cdef class DPosition1:
    cdef shared_ptr[_DPosition1] inst
cdef class DPosition2:
    cdef shared_ptr[_DPosition2] inst


cdef class DPosition1:

    def __dealloc__(self):
        self.inst.reset()

    def __init__(self, *args , **kwargs):
        if not args:
             self._init_0(*args)
        elif (len(args)==1):
             self._init_1(*args)
        else:
             raise Exception('can not handle type of %s' % (args,))
 
    def _init_0(self):
        self.inst = shared_ptr[_DPosition1](new _DPosition1())

    def _init_1(self, double a):
        self.inst = shared_ptr[_DPosition1](new _DPosition1(a))

    def __getitem__(self, ix):
        if ix != 0:
            raise IndexError("invalid index %d" % ix)
        return deref(self.inst.get())[0]

    def __repr__(self):
        """
        __repr__(self: DPosition1) -> str
        
        Return a string representation of the DPosition1 object.

        Returns key properties in a readable format:
        DPosition1(x=100.0)
        """
        cdef double x = deref(self.inst.get())[0]
        return f"DPosition1(x={x:.4f})"


cdef class DPosition2:

    def __dealloc__(self):
        self.inst.reset()


    def __init__(self, *args , **kwargs):
        if not args:
             self._init_0(*args)
        elif (len(args)==2):
             self._init_1(*args)
        else:
             raise Exception('can not handle type of %s' % (args,))
 
    def _init_0(self):
        self.inst = shared_ptr[_DPosition2](new _DPosition2())

    def _init_1(self, double a, double b):
        self.inst = shared_ptr[_DPosition2](new _DPosition2(a, b))

    def __getitem__(self, ix):
        if ix != 0 and ix != 1:
            raise IndexError("invalid index %d" % ix)
        return deref(self.inst.get())[ix]

    def __repr__(self):
        """
        __repr__(self: DPosition2) -> str
        
        Return a string representation of the DPosition2 object.

        Returns key properties in a readable format:
        DPosition2(x=100.0, y=200.0)
        """
        cdef double x = deref(self.inst.get())[0]
        cdef double y = deref(self.inst.get())[1]
        return f"DPosition2(x={x:.4f}, y={y:.4f})"
