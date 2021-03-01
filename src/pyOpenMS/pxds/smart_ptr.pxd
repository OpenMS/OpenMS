
#cdef extern from "boost/smart_ptr/shared_ptr.hpp" namespace "boost":
cdef extern from "<memory>" namespace "std":

    cdef cppclass shared_ptr[T]:
        shared_ptr()
        shared_ptr(T*)
        void reset()
        T* get()
        int unique()
        int use_count()
