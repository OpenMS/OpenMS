from Matrix cimport *
cimport numpy as np
from numpy.lib.stride_tricks import as_strided


# continue with extra code if needed
            

    def get_matrix_as_view(self):
        """Cython signature: numpy_matrix get_matrix_as_view()
        """

        cdef _Matrix[double] * mat_ = self.inst.get()
        cdef unsigned int rows = mat_.rows()
        cdef unsigned int cols = mat_.cols()
        cdef double* data = mat_.data()
        cdef double[:,:] mem_view = <double[:rows,:cols]>data
        dtype = 'double'
        cdef int itemsize = np.dtype(dtype).itemsize
        cdef unsigned int row_stride, col_stride
        if mat_.rowMajor():
            row_stride = mat_.outerStride() if mat_.outerStride() > 0 else cols
            col_stride = mat_.innerStride() if mat_.innerStride() > 0 else 1
        else:
            row_stride = mat_.innerStride() if mat_.innerStride() > 0 else 1
            col_stride = mat_.outerStride() if mat_.outerStride() > 0 else rows

        return np.lib.stride_tricks.as_strided(np.asarray(mem_view, dtype=dtype, order="F"), strides=[row_stride*itemsize, col_stride*itemsize])


    def get_matrix(self):
        """Cython signature: numpy_matrix get_matrix()
        """

        cdef _Matrix[double] * mat_ = self.inst.get()
        cdef unsigned int rows = mat_.rows()
        cdef unsigned int cols = mat_.cols()
        cdef double* data = mat_.data()
        cdef double[:,:] mem_view = <double[:rows,:cols]>data
        dtype = 'double'
        cdef int itemsize = np.dtype(dtype).itemsize
        cdef unsigned int row_stride, col_stride
        if mat_.rowMajor():
            row_stride = mat_.outerStride() if mat_.outerStride() > 0 else cols
            col_stride = mat_.innerStride() if mat_.innerStride() > 0 else 1
        else:
            row_stride = mat_.innerStride() if mat_.innerStride() > 0 else 1
            col_stride = mat_.outerStride() if mat_.outerStride() > 0 else rows

        return np.copy(np.lib.stride_tricks.as_strided(np.asarray(mem_view, dtype=dtype, order="F"), strides=[row_stride*itemsize, col_stride*itemsize]))

#    def set_matrix(self, np.ndarray[double, ndim=2, mode="c"] data not None):
#        """Cython signature: numpy_matrix set_matrix()
#        """

#        cdef _Matrix[double] * mat_ = self.inst.get()

#        cdef unsigned int rows = data.shape[0]
#        cdef unsigned int cols = data.shape[1]
#        mat_.resize(rows, cols)

#        cdef int i = 0
#        cdef int j = 0
#        for i in range(int(rows)):
#            for j in range(int(cols)):
#                mat_.setValue(i,j,data[i][j])

