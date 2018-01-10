

    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_RNaseDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_RNaseDB](_getInstance_RNaseDB())
