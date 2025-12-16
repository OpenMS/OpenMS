from ConvexHull2D cimport ConvexHull2D as _ConvexHull2D
from DPosition cimport DPosition2 as _DPosition2
from DBoundingBox cimport DBoundingBox2 as _DBoundingBox2


    def enclosesXY(self, float x, float y):
        """
        enclosesXY(self: ConvexHull2D, x: float, y: float) -> int
        
        Check if a point (x, y) is enclosed by the convex hull.
        
        Parameters:
        x (float)
        y (float)
        
        Returns:
        int
        """
        cdef _DPosition2 pos
        pos[0] = x
        pos[1] = y
        return self.inst.get().encloses(pos)


    def getHullPointsNPY(self):
        """
        getHullPointsNPY(self: ConvexHull2D) -> np.ndarray
        
        Get the hull points as a numpy array.
        
        Returns:
        result (np.ndarray[np.float32_t, ndim=2])
        """
        cdef libcpp_vector[_DPosition2] points = self.inst.get().getHullPoints()
        cdef np.ndarray[np.float32_t, ndim=2] result
        cdef n = points.size()
        result = np.zeros( [n,2], dtype=np.float32)
        cdef libcpp_vector[_DPosition2].iterator it = points.begin()
        cdef int i = 0
        while it != points.end():
            result[i,0] = deref(it)[0]
            result[i,1] = deref(it)[1]
            inc(it)
            i += 1
        return result

    def setHullPointsNPY(self, np.ndarray[np.float32_t, ndim=2] points):
        """
        setHullPointsNPY(self: ConvexHull2D, points: np.ndarray) -> None
        
        Set the hull points from a numpy array.
        
        Parameters:
        points (np.ndarray[np.float32_t, ndim=2])
        """
        cdef _ConvexHull2D * hull = self.inst.get()
        cdef int N = points.shape[0]
        cdef int i
        cdef libcpp_vector[_DPosition2] vec
        cdef _DPosition2 p
        for i in range(N):
            p[0] = points[i,0]
            p[1] = points[i,1]
            vec.push_back(p)
        self.inst.get().setHullPoints(vec)

    def getBoundingBox2D(self):
        """
        getBoundingBox2D(self: ConvexHull2D) -> Tuple[Tuple[float, float], Tuple[float, float]]
        
        Get the bounding box of the convex hull.
        
        Returns:
        ((double,double),(double,double))
        """
        cdef _DBoundingBox2 box = self.inst.get().getBoundingBox()
        cdef _DPosition2 minp = box.minPosition()
        cdef _DPosition2 maxp = box.maxPosition()
        return (minp[0], minp[1]), (maxp[0], maxp[1])

    def addPointXY(self, x, y):
        """
        addPointXY(self: ConvexHull2D, x: float, y: float) -> None
        
        Add a point to the convex hull.
        
        Parameters:
        x (double)
        y (double)
        """
        cdef _DPosition2 p
        p[0] = x
        p[1] = y
        self.inst.get().addPoint(p)

    def addPointsNPY(self, np.ndarray[np.float32_t, ndim=2] points):
        """
        addPointsNPY(self: ConvexHull2D, points: np.ndarray) -> None
        
        Add multiple points to the convex hull.
        
        Parameters:
        points (np.ndarray[np.float32_t, ndim=2])
        """
        cdef _ConvexHull2D * hull = self.inst.get()
        cdef int N = points.shape[0]
        cdef int i
        cdef libcpp_vector[_DPosition2] vec
        cdef _DPosition2 p
        for i in range(N):
            p[0] = points[i,0]
            p[1] = points[i,1]
            vec.push_back(p)
        self.inst.get().addPoints(vec)

    def __repr__(self):
        """
        __repr__(self: ConvexHull2D) -> str
        
        Return a string representation of the ConvexHull2D object.

        Returns key properties in a readable format:
        ConvexHull2D(num_points=10, bbox=((100.0, 400.0), (110.0, 410.0)))
        """
        cdef libcpp_vector[_DPosition2] points = self.inst.get().getHullPoints()
        cdef int num_points = points.size()
        
        parts = []
        parts.append(f"num_points={num_points}")
        
        # Get bounding box if we have points (wrap in try/except for safety)
        if num_points > 0:
            try:
                bbox = self.getBoundingBox2D()
                parts.append(f"bbox={bbox}")
            except:
                pass  # Skip bounding box if it can't be computed
        
        return f"ConvexHull2D({', '.join(parts)})"
