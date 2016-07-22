# This will go into the header
from ElementDB cimport ElementDB as _ElementDB
from Map cimport Map as _Map
cdef class ElementDBWrapper:
    # A small utility class holding a ptr and implementing get()
    cdef const _ElementDB* wrapped
    cdef setptr(self, const _ElementDB* wrapped): self.wrapped = wrapped
    cdef const _ElementDB* get(self) except *: return self.wrapped


    # This will go into the class
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    cdef ElementDBWrapper inst

    def __init__(self):
      self.inst = ElementDBWrapper()
      self.inst.setptr(_getInstance_ElementDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass

    # def getAtomicNumbers(self):
    #     _r = self.inst.get().getAtomicNumbers()
    #     py_result = dict()
    #     cdef _Map[unsigned int, _Element *].iterator it__r = _r.begin()
    #     cdef Element item_py_result
    #     while it__r != _r.end():
    #        item_py_result = Element.__new__(Element)
    #        # item_py_result.inst = shared_ptr[_Element *](new _Element *((deref(it__r)).second))
    #        item_py_result.inst = shared_ptr[_Element](new _Element(deref(deref(it__r).second)))
    #        py_result[<unsigned int>(deref(it__r).first)] = item_py_result
    #        inc(it__r)
    #     return py_result

