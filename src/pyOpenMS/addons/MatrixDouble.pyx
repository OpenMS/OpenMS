from Matrix cimport *


# continue with extra code if needed

    def get_matrix(self):
        """Cython signature: numpy_matrix get_matrix()
        """

        cdef _Matrix[double] * mat_ = self.inst.get()

        cdef unsigned int rows = mat_.rows()
        cdef unsigned int cols = mat_.cols()

        cdef libcpp_vector[double] tmp_vec;
        tmp_vec = mat_.asVector();

        xarr = np.asarray(tmp_vec)
        xarr = xarr.reshape(rows, cols)
        return xarr

    def get_matrix_as_view(self):
        """Cython signature: numpy_matrix get_matrix()
        """

        cdef _Matrix[double] * mat_ = self.inst.get()

        cdef unsigned int rows = mat_.rows()
        cdef unsigned int cols = mat_.cols()
        cdef unsigned int n = rows * cols
        cdef np.ndarray[double, ndim=2] data
        data = np.zeros( (rows,cols), dtype=np.float64)

        cdef libcpp_vector[double] * vec_ptr = <libcpp_vector[double]*> mat_
        cdef double * raw_ptr =  address(deref(vec_ptr)[0]) 

        ## # We use a memory view to get the data from the raw data
        ## # See https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html 
        ## # See https://stackoverflow.com/questions/43021574/cast-c-array-into-numpy-array-cython-typed-memoryview-in-cython-code
        cdef double[:] vec_view = <double[:n]>raw_ptr # cast to memoryview, refer to the underlying buffer without copy
        xarr = np.asarray(vec_view) # numpy array refer to the underlying buffer without copy
        xarr = xarr.reshape(rows, cols)
        return xarr

    def set_matrix(self, np.ndarray[double, ndim=2, mode="c"] data not None):
        """Cython signature: numpy_matrix set_matrix()
        """

        cdef _Matrix[double] * mat_ = self.inst.get()

        cdef unsigned int rows = data.shape[0]
        cdef unsigned int cols = data.shape[1]
        mat_.resize(rows, cols, 0)

        cdef int i = 0
        cdef int j = 0
        for i in range(int(rows)):
            for j in range(int(cols)):
                mat_.setValue(i,j,data[i][j])



