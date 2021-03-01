
cdef extern from "<memory>" namespace "std":

    cdef cppclass shared_ptr[T]:
        shared_ptr()
        shared_ptr(T*)
        void reset()
        T* get()
        int unique()
        int use_count()
