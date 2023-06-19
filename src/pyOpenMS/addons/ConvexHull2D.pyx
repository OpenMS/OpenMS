from ConvexHull2D cimport ConvexHull2D as _ConvexHull2D
from DPosition cimport DPosition2 as _DPosition2
from DBoundingBox cimport DBoundingBox2 as _DBoundingBox2


    def enclosesXY(self, float x, float y):
        """
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
        Returns:
        ((double,double),(double,double))
        """
        cdef _DBoundingBox2 box = self.inst.get().getBoundingBox()
        cdef _DPosition2 minp = box.minPosition()
        cdef _DPosition2 maxp = box.maxPosition()
        return (minp[0], minp[1]), (maxp[0], maxp[1])

    def addPointXY(self, x, y):
        """
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

