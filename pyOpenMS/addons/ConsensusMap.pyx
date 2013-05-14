

    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def getFileDescriptions(self):

        # FileDescriptions is a type alias for Map<..> which can not be
        # handled by autowrap automatically. So we have to provide a manual
        # converter here:
        #
        # (as far as I remember the wrapper works on linux, but msvc complains
        # a lot about the generated code... uwe schmitt)

        cdef FileDescriptions _r = self.inst.get().getFileDescriptions()
        py_result = dict()
        cdef FileDescriptions_iterator it__r = _r.begin()
        cdef FileDescription item_py_result
        while it__r != _r.end():
           item_py_result = FileDescription.__new__(FileDescription)
           item_py_result.inst = shared_ptr[_FileDescription](new _FileDescription((deref(it__r)).second))
           py_result[<unsigned long int>(deref(it__r).first)] = item_py_result
           inc(it__r)
        return py_result

    def setFileDescriptions(self, dict in_0 ):
        assert isinstance(in_0, dict) and all(isinstance(k, (int, long)) for k in in_0.keys()) and all(isinstance(v, FileDescription) for v in in_0.values()), 'arg in_0 wrong type'
        cdef FileDescriptions v0
        for key, value in in_0.items():
           v0[<unsigned long int> key] = deref((<FileDescription>value).inst.get())

        # we have to utilize AutowrapRefHolder here, because cython does not
        # like
        #
        #     self.inst.get().getFileDescriptions() = v0
        #
        # and if you choose
        #
        #     cdef FileDescriptions & ref = self.inst.get().getFileDescriptions()
        #     ref = v0
        #
        # the c++ compiler will complain, as cython generates invalid code in
        # this case: Cython creates two statements from the first line:
        #
        #     FileDescription & ref;   // invalid: ref not initailized
        #     ref = self.inst.get().getFileDescriptions()
        #
        # wrapping this ref into a class holding the ref and using pointers
        # works. the following code results in c++ code simialar to
        #
        #    AutowrapRefHolder<FileDescription> * refholder;
        #    refholder = new AutowrapRefHolder<FileDescription>(self.inst...)
        #    refholder.assign(v0)

        cdef AutowrapRefHolder[FileDescriptions] * refholder
        refholder = new AutowrapRefHolder[FileDescriptions](self.inst.get().getFileDescriptions())
        refholder.assign(v0)

        cdef replace_in_0 = dict()
        cdef FileDescriptions_iterator it_in_0 = v0.begin()
        cdef FileDescription item_in_0
        while it_in_0 != v0.end():
           item_in_0 = FileDescription.__new__(FileDescription)
           item_in_0.inst = shared_ptr[_FileDescription](new _FileDescription((deref(it_in_0)).second))
           replace_in_0[<unsigned long int> deref(it_in_0).first] = item_in_0
           inc(it_in_0)
        in_0.clear()
        in_0.update(replace_in_0)
