

    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_RNaseDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_RNaseDB](_getInstance_RNaseDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass
