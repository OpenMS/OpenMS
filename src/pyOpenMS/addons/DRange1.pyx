


    def _init_2(self, double min_val, double max_val):
        """
        _init_2(self, min_val: float, max_val: float) -> None

        Initialize DRange1 with min and max values.
        """
        self.inst = shared_ptr[_DRange1](new _DRange1())
        self.inst.get().setMinX(min_val)
        self.inst.get().setMaxX(max_val)

    def __init__(self, *args, **kwargs):
        """
        __init__(self, *args) -> None

        Overloaded constructors:
        - DRange1() - creates an empty range
        - DRange1(other: DRange1) - copy constructor
        - DRange1(min: float, max: float) - creates a range with given bounds
        """
        if not args:
            self._init_0(*args)
        elif (len(args) == 1) and (isinstance(args[0], DRange1)):
            self._init_1(*args)
        elif (len(args) == 2) and (isinstance(args[0], (int, float))) and (isinstance(args[1], (int, float))):
            self._init_2(float(args[0]), float(args[1]))
        else:
            raise Exception('can not handle type of %s' % (args,))

    def getMin(self):
        """
        getMin(self) -> float

        Returns the minimum value of this 1D range.

        Returns:
            float: The minimum x coordinate
        """
        return self.minX()

    def getMax(self):
        """
        getMax(self) -> float

        Returns the maximum value of this 1D range.

        Returns:
            float: The maximum x coordinate
        """
        return self.maxX()

    def encloses(self, double position):
        """
        encloses(self, position: float) -> bool

        Check if a position is within this range.
        The range is [min, max), so min is included but max is excluded.

        Args:
            position: The position to check

        Returns:
            bool: True if the position is within the range
        """
        cdef double min_val = self.minX()
        cdef double max_val = self.maxX()
        return position >= min_val and position < max_val

    def __repr__(self):
        """
        __repr__(self) -> str

        Return a string representation of the DRange1 object.
        """
        cdef double min_val = self.minX()
        cdef double max_val = self.maxX()
        return f"DRange1([{min_val}, {max_val}))"

    def __str__(self):
        """
        __str__(self) -> str

        Return a string representation of the DRange1 object.
        """
        return self.__repr__()
