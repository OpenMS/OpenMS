


    def setModel(self, PILISModel model):
        self.inst.get().setModel( model.inst.get() )

    def getIdentification(self, dict candidates, PeptideIdentification id_, RichMSSpectrum spectrum):
        cdef libcpp_map[_String, unsigned int] c_dict
        for k,v in candidates.iteritems():
            c_dict[ _String(<char *>k) ] = v
        
        self.inst.get().getIdentification(c_dict, deref(id_.inst.get()), deref(spectrum.inst.get()) )
