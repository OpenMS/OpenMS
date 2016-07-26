

    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    cdef AutowrapPtrHolder[_EnzymesDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_EnzymesDB](_getInstance_EnzymesDB())

