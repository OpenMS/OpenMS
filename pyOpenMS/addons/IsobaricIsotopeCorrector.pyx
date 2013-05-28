


    def __init__ItraqEightPlexQuantitationMethod(self, ItraqEightPlexQuantitationMethod this_quant_method):
        self.initialize_from_ptr( this_quant_method.inst.get() )

    def __init__ItraqFourPlexQuantitationMethod(self, ItraqFourPlexQuantitationMethod this_quant_method):
        self.initialize_from_ptr( this_quant_method.inst.get() )

    def __init__TMTSixPlexQuantitationMethod(self, TMTSixPlexQuantitationMethod this_quant_method):
        self.initialize_from_ptr( this_quant_method.inst.get() )

    cdef initialize_from_ptr(self, _IsobaricQuantitationMethod * this_quant_method):
        self.inst = shared_ptr[_IsobaricIsotopeCorrector](new _IsobaricIsotopeCorrector(( this_quant_method )))

    def __init__from_self(self, IsobaricIsotopeCorrector in_0 ):
        assert isinstance(in_0, IsobaricIsotopeCorrector), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsobaricIsotopeCorrector](new _IsobaricIsotopeCorrector((deref(in_0.inst.get()))))

    def __init__(self, input_arg):
        if isinstance(input_arg, IsobaricIsotopeCorrector):
            self.__init__from_self(input_arg)
        elif isinstance(input_arg, ItraqEightPlexQuantitationMethod):
            self.__init__ItraqEightPlexQuantitationMethod(input_arg)
        elif isinstance(input_arg, ItraqFourPlexQuantitationMethod):
            self.__init__ItraqFourPlexQuantitationMethod(input_arg)
        elif isinstance(input_arg, TMTSixPlexQuantitationMethod):
            self.__init__TMTSixPlexQuantitationMethod(input_arg)
        else:
            raise Exception("Cannot create with input type %s" % type(input_arg))
