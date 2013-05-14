

    def fromExperiment(self, MSExperiment in_0 ):
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'

        self.inst = shared_ptr[_MSExperiment[_Peak1D,_ChromatogramPeak]](new _MSExperiment[_Peak1D,_ChromatogramPeak]( <_MSExperiment[_Peak1D,_ChromatogramPeak] &>deref(in_0.inst.get())  ))

