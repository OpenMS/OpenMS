

    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # Provides the same interface as expected by autowrap (get/assign) but
    # simply holds the ptr and does not do any memory management.
    # see autowrap/data_files/autowrap/README.md for implementation
    # cdef AutowrapPtrHolder[_CrossLinksDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_CrossLinksDB](_getInstance_CrossLinksDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass


