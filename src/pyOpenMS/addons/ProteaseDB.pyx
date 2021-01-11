

    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_ProteaseDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_ProteaseDB](_getInstance_ProteaseDB())
