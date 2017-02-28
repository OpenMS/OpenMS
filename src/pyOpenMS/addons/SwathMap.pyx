cdef extern from *:
    cdef T dynamic_cast[T](void *) except +
cdef extern from "boost/smart_ptr/shared_ptr.hpp" namespace "boost":
    cdef shared_ptr[T] dynamic_pointer_cast[T,U](shared_ptr[U]& r) except +
from SpectrumAccessOpenMS cimport SpectrumAccessOpenMS as _SpectrumAccessOpenMS
from SpectrumAccessOpenMSCached cimport SpectrumAccessOpenMSCached as _SpectrumAccessOpenMSCached
from SpectrumAccessOpenMSInMemory cimport SpectrumAccessOpenMSInMemory as _SpectrumAccessOpenMSInMemory
from SpectrumAccessQuadMZTransforming cimport SpectrumAccessQuadMZTransforming as _SpectrumAccessQuadMZTransforming
ctypedef _SpectrumAccessOpenMS* _SpectrumAccessOpenMSPtr
ctypedef _SpectrumAccessOpenMSCached* _SpectrumAccessOpenMSCachedPtr
ctypedef _SpectrumAccessOpenMSInMemory* _SpectrumAccessOpenMSInMemoryPtr
ctypedef _SpectrumAccessQuadMZTransforming * _SpectrumAccessQuadMZTransformingPtr


    def getSpectrumPtr(self):
        _r = self.inst.get().sptr

        if (_r.get() == NULL):
          return None

        cdef _SpectrumAccessOpenMS * ptr_sa = dynamic_cast[ _SpectrumAccessOpenMSPtr ](_r.get() )
        cdef _SpectrumAccessOpenMSInMemory * ptr_inmem = dynamic_cast[ _SpectrumAccessOpenMSInMemoryPtr ](_r.get() )
        cdef _SpectrumAccessOpenMSCached * ptr_cached = dynamic_cast[ _SpectrumAccessOpenMSCachedPtr ](_r.get() )
        cdef _SpectrumAccessQuadMZTransforming * ptr_quad = dynamic_cast[ _SpectrumAccessQuadMZTransformingPtr ](_r.get() )

        if (ptr_sa != NULL):
          res_sa = SpectrumAccessOpenMS(__createUnsafeObject__=True)
          res_sa.inst = dynamic_pointer_cast[_SpectrumAccessOpenMS, _ISpectrumAccess](_r)
          return res_sa
        elif (ptr_inmem != NULL):
          res_inmem = SpectrumAccessOpenMSInMemory(__createUnsafeObject__=True)
          res_inmem.inst = dynamic_pointer_cast[_SpectrumAccessOpenMSInMemory, _ISpectrumAccess](_r)
          return res_inmem
        elif (ptr_cached != NULL):
          res_cached = SpectrumAccessOpenMSCached(__createUnsafeObject__=True)
          res_cached.inst = dynamic_pointer_cast[_SpectrumAccessOpenMSCached, _ISpectrumAccess](_r)
          return res_cached
        elif (ptr_quad != NULL):
          res_quad = SpectrumAccessQuadMZTransforming(__createUnsafeObject__=True)
          res_quad.inst = dynamic_pointer_cast[_SpectrumAccessQuadMZTransforming, _ISpectrumAccess](_r)
          return res_quad
        else:
          raise Exception("Did not find suitable conversion to Python object")

    def setSpectrumPtr(self, arg):
        cdef SpectrumAccessOpenMS arg_sa 
        cdef SpectrumAccessOpenMSCached arg_cached 
        cdef SpectrumAccessOpenMSInMemory arg_inmem 
        cdef SpectrumAccessQuadMZTransforming arg_quad 
        if isinstance(arg, SpectrumAccessOpenMS):
            arg_sa = arg
            self.inst.get().sptr = dynamic_pointer_cast[_ISpectrumAccess,_SpectrumAccessOpenMS](arg_sa.inst)
        elif isinstance(arg, SpectrumAccessOpenMSCached):
            arg_cached = arg
            self.inst.get().sptr = dynamic_pointer_cast[_ISpectrumAccess,_SpectrumAccessOpenMSCached](arg_cached.inst)
        elif isinstance(arg, SpectrumAccessOpenMSInMemory):
            arg_inmem = arg
            self.inst.get().sptr = dynamic_pointer_cast[_ISpectrumAccess,_SpectrumAccessOpenMSInMemory](arg_inmem.inst)
        elif isinstance(arg, SpectrumAccessQuadMZTransforming):
            arg_quad = arg
            self.inst.get().sptr = dynamic_pointer_cast[_ISpectrumAccess,_SpectrumAccessQuadMZTransforming](arg_quad.inst)
        else:
          raise Exception("Need to provide suitable ISpectrumAccess-derived child class")

