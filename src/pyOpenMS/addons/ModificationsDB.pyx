# This will go into the header
from ModificationsDB cimport ModificationsDB as _ModificationsDB
from Map cimport Map as _Map
cdef class ModificationsDBWrapper:
    # A small utility class holding a ptr and implementing get()
    cdef _ModificationsDB* wrapped
    cdef setptr(self, _ModificationsDB* wrapped): self.wrapped = wrapped
    cdef _ModificationsDB* get(self) except *: return self.wrapped


    # This will go into the class
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    cdef ModificationsDBWrapper inst

    def __init__(self):
      self.inst = ModificationsDBWrapper()
      self.inst.setptr(_getInstance_ModificationsDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass


