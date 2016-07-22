# This will go into the header
from EnzymesDB cimport EnzymesDB as _EnzymesDB
cdef class EnzymesDBWrapper:
    # A small utility class holding a ptr and implementing get()
    cdef _EnzymesDB* wrapped
    cdef setptr(self, _EnzymesDB* wrapped): self.wrapped = wrapped
    cdef _EnzymesDB* get(self) except *: return self.wrapped


    # This will go into the class
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    cdef EnzymesDBWrapper inst

    def __init__(self):
      self.inst = EnzymesDBWrapper()
      self.inst.setptr(_getInstance_EnzymesDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass

