

    def __int__(self):
        # <long> calls c++ conversion streampos -> c++ long value
        # int(..) converts this value to a Python long.
        # value. else this method may return an int or long, depending on the magnitude of
        # the value !
        return int(<long>(deref(self.inst.get())))
