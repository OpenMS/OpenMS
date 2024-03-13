from Matrix cimport *
cimport numpy as np

# continue with extra code if needed

    def get_matrix_as_view(self):
        """Cython signature: numpy_matrix get_matrix_as_view()
        """

        cdef _Matrix[double] * mat_ = self.inst.get()

        cdef unsigned int rows = mat_.rows()
        cdef unsigned int cols = mat_.cols()
        cdef double* data = mat_.data()
        # Create a NumPy array from the buffer, specifying not to copy the data.
        return np.array(<double[:rows * cols:1]> data, copy=False).reshape(rows, cols)


    def get_matrix(self):
        """Cython signature: numpy_matrix get_matrix()
        """

        cdef _Matrix[double] * mat_ = self.inst.get()
        cdef unsigned int rows = mat_.rows()
        cdef unsigned int cols = mat_.cols()
        cdef double* data = mat_.data()
        # Create a NumPy array from the buffer, specifying not to copy the data.
        return np.array(<double[:rows * cols:1]> data, copy=True).reshape(rows, cols)


#    def set_matrix(self, np.ndarray[double, ndim=2, mode="c"] data not None):
#        """Cython signature: numpy_matrix set_matrix()
#        """
#
#        cdef _Matrix[double] * mat_ = self.inst.get()
#
#        cdef unsigned int rows = data.shape[0]
#        cdef unsigned int cols = data.shape[1]
#        mat_.resize(rows, cols, 0)
#
#        cdef int i = 0
#        cdef int j = 0
#        for i in range(int(rows)):
#            for j in range(int(cols)):
#                mat_.setValue(i,j,data[i][j])



