from Matrix cimport *
cimport numpy as np
from numpy.lib.stride_tricks import as_strided


# continue with extra code if needed
            

    
    @staticmethod
    def fromNdArray(np.ndarray[double, ndim=2] data not None):
        """Creates a new Matrix from a numpy ndarray."""
        cdef MatrixDouble mat = MatrixDouble()
        mat.set_matrix(data)
        return mat

    def get_matrix_as_view(self):
        """get_matrix(self) -> np.ndarray[double, ndim=2]

        Returns a view on the underlying Matrix as a 2D numpy ndarray.
        .. caution::
           Future changes to the Matrix will affect the ndarray and vice versa.
           Make sure that the Matrix does not go out of scope before the last use
           of your ndarray.
        """

        cdef _Matrix[double] * mat_ = self.inst.get()
        cdef unsigned int rows = mat_.rows()
        cdef unsigned int cols = mat_.cols()
        cdef double* data = mat_.data()
        cdef double[:,:] mem_view = <double[:rows,:cols]>data
        dtype = 'double'
        cdef int itemsize = np.dtype(dtype).itemsize
        cdef unsigned int row_stride, col_stride
        o = 'F'
        if mat_.rowMajor():
            row_stride = mat_.outerStride() if mat_.outerStride() > 0 else cols
            col_stride = mat_.innerStride() if mat_.innerStride() > 0 else 1
            o = 'F'
        else:
            row_stride = mat_.innerStride() if mat_.innerStride() > 0 else 1
            col_stride = mat_.outerStride() if mat_.outerStride() > 0 else rows
            o = 'C'

        return np.lib.stride_tricks.as_strided(np.asarray(mem_view, dtype=dtype, order=o), strides=[row_stride*itemsize, col_stride*itemsize])


    def get_matrix(self):
        """get_matrix(self) -> np.ndarray[double, ndim=2]

        Returns a copy of the underlying Matrix as a 2D numpy ndarray.
        """

        return np.copy(self.get_matrix_as_view())


    def set_matrix(self, np.ndarray[double, ndim=2] data not None):
        """set_matrix(self, data: np.ndarray[double, ndim=2]) -> None

        Copies the values from the numpy ndarray into the Matrix.
        """

        cdef _Matrix[double] * mat_ = self.inst.get()

        cdef unsigned int rows = data.shape[0]
        cdef unsigned int cols = data.shape[1]
        mat_.resize(rows, cols)

        cdef int i = 0
        cdef int j = 0
        for i in range(int(rows)):
            for j in range(int(cols)):
                mat_.setValue(i, j, data[i][j])

