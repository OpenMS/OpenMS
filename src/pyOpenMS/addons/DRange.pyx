# This goes into the PXD file
from DRange cimport DRange1 as _DRange1
from DPosition cimport DPosition1 as _DPosition1
cdef class DRange1:
    cdef shared_ptr[_DRange1] inst


cdef class DRange1:
    """
    A 1-dimensional range (interval) for filtering by a single dimension such as RT, m/z, or intensity.

    This class represents a half-open interval [min, max) in one dimension.
    It is commonly used with PeakFileOptions and FeatureFileOptions to filter
    data by retention time, m/z, or intensity ranges.

    Example usage:
        >>> rt_range = DRange1(100.0, 200.0)  # RT range from 100 to 200
        >>> options = PeakFileOptions()
        >>> options.setRTRange(rt_range)
    """

    def __dealloc__(self):
        self.inst.reset()

    def __init__(self, *args, **kwargs):
        """
        Create a DRange1 object.

        Args:
            No arguments: Creates an empty range
            Two float arguments (min_val, max_val): Creates a range [min_val, max_val)

        Examples:
            >>> r = DRange1()           # empty range
            >>> r = DRange1(10.0, 20.0) # range from 10 to 20
        """
        if not args:
            self._init_0(*args)
        elif len(args) == 2:
            self._init_2(*args)
        else:
            raise Exception('DRange1 takes either 0 or 2 arguments, got %d: %s' % (len(args), args))

    def _init_0(self):
        self.inst = shared_ptr[_DRange1](new _DRange1())

    def _init_2(self, double min_val, double max_val):
        cdef _DPosition1 min_pos = _DPosition1(min_val)
        cdef _DPosition1 max_pos = _DPosition1(max_val)
        self.inst = shared_ptr[_DRange1](new _DRange1(min_pos, max_pos))

    @property
    def minPosition(self):
        """Returns the minimum value of the range."""
        cdef _DPosition1 pos = deref(self.inst.get()).minPosition()
        return pos[0]

    @property
    def maxPosition(self):
        """Returns the maximum value of the range."""
        cdef _DPosition1 pos = deref(self.inst.get()).maxPosition()
        return pos[0]

    def setMin(self, double val):
        """Set the minimum value of the range."""
        cdef _DPosition1 pos = _DPosition1(val)
        deref(self.inst.get()).setMin(pos)

    def setMax(self, double val):
        """Set the maximum value of the range."""
        cdef _DPosition1 pos = _DPosition1(val)
        deref(self.inst.get()).setMax(pos)

    def setMinMax(self, double min_val, double max_val):
        """Set both minimum and maximum values of the range."""
        cdef _DPosition1 min_pos = _DPosition1(min_val)
        cdef _DPosition1 max_pos = _DPosition1(max_val)
        deref(self.inst.get()).setMinMax(min_pos, max_pos)

    def __repr__(self):
        return "DRange1([%s, %s))" % (self.minPosition, self.maxPosition)

    def __str__(self):
        return self.__repr__()

