



    def getData(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        retval = np.asarray(arr.copy())
        return retval

    def getData_mv(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        return arr
