# This will go into the header
from ResidueDB cimport ResidueDB as _ResidueDB
from Map cimport Map as _Map
cdef class ResidueDBWrapper:
    # A small utility class holding a ptr and implementing get()
    cdef _ResidueDB* wrapped
    cdef setptr(self, _ResidueDB* wrapped): self.wrapped = wrapped
    cdef _ResidueDB* get(self) except *: return self.wrapped


    # This will go into the class
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    cdef ResidueDBWrapper inst

    def __init__(self):
      self.inst = ResidueDBWrapper()
      self.inst.setptr(_getInstance_ResidueDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass

