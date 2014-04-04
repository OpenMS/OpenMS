

    def __long__(self):
        # <long> calls c++ conversion streampos -> c++ long value
        # long(..) converts this value to a Python long.
        # value. else this method may return an int or long, depending on the magnitude of
        # the value !
        return long(<long>(deref(self.inst.get())))
