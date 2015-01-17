


    ## def __extractChromatograms_OpenMS(self, SpectrumAccessOpenMS input,
    ##         libcpp_vector[shared_ptr[_OSChromatogram]] & v1,
    ##         libcpp_vector[_ExtractionCoordinates] & v2,
    ##         double mz_extraction_window, bool ppm, _String filter):

    ##     self._extractChromatograms_inner( this_quant_method.inst.get() )

    ## cdef _extractChromatograms_inner(self,
    ##         shared_ptr[ _ISpectrumAccess ] & input,
    ##         libcpp_vector[shared_ptr[_OSChromatogram]] & v1,
    ##         libcpp_vector[_ExtractionCoordinates] & v2,
    ##         double mz_extraction_window, bool ppm, _String filter):

    ##     #self.inst = shared_ptr[_IsobaricChannelExtractor](new _IsobaricChannelExtractor(( this_quant_method )))

    ##     self.inst.get().extractChromatograms(input, v1, deref(v2), (<double>mz_extraction_window), (<bool>ppm), (_String(<char *>filter)))



    #  ### This doesnt work,
    #  ###  cannot initialize a member subobject of type 'OpenSwath::ISpectrumAccess *' with an lvalue of type '_object *'
    #  ###  
    #  def extractChromatograms(self, input , list output , list extraction_coordinates , double mz_extraction_window ,  ppm , bytes filter ):
    #      # assert isinstance(input, SpectrumAccessOpenMS), 'arg input wrong type'
    #      assert isinstance(output, list) and all(isinstance(elemt_rec, OSChromatogram) for elemt_rec in output), 'arg output wrong type'
    #      assert isinstance(extraction_coordinates, list) and all(isinstance(elemt_rec, ExtractionCoordinates) for elemt_rec in extraction_coordinates), 'arg extraction_coordinates wrong type'
    #      assert isinstance(mz_extraction_window, float), 'arg mz_extraction_window wrong type'
    #      assert isinstance(ppm, (int, long)), 'arg ppm wrong type'
    #      assert isinstance(filter, bytes), 'arg filter wrong type'

    #      # cdef shared_ptr[_SpectrumAccessOpenMS] input_input = input.inst
    #      cdef libcpp_vector[shared_ptr[_OSChromatogram]] v1
    #      cdef OSChromatogram output_rec
    #      for output_rec in output:
    #          v1.push_back(output_rec.inst)
    #      cdef libcpp_vector[_ExtractionCoordinates] * v2 = new libcpp_vector[_ExtractionCoordinates]()
    #      cdef ExtractionCoordinates item2
    #      for item2 in extraction_coordinates:
    #          v2.push_back(deref(item2.inst.get()))
    #  
    #      # Call 
    #      self.inst.get().extractChromatograms((<shared_ptr[_ISpectrumAccess]>input.inst), v1, deref(v2), (<double>mz_extraction_window), (<bool>ppm), (_String(<char *>filter)))

    #      del v2

    #      # gather results
    #      replace = list()
    #      cdef libcpp_vector[shared_ptr[_OSChromatogram]].iterator it_output = v1.begin()
    #      while it_output != v1.end():
    #         output_rec = OSChromatogram.__new__(OSChromatogram)
    #         output_rec.inst = deref(it_output)
    #         replace.append(output_rec)
    #         inc(it_output)
    #      # replace old vector with new contents
    #      output[:] = []
    #      for output_rec_b in replace:
    #          output.append(output_rec_b)

    def extractChromatograms(self, input , list output , list extraction_coordinates , double mz_extraction_window ,  ppm , bytes filter ):
        if isinstance(input, SpectrumAccessOpenMS):
            self._extractChromatograms_OpenMS(input, output, extraction_coordinates, mz_extraction_window, ppm, filter)
        elif isinstance(input, SpectrumAccessOpenMSCached):
            self._extractChromatograms_OpenMSCached(input, output, extraction_coordinates, mz_extraction_window, ppm, filter)

    def _extractChromatograms_OpenMS(self, SpectrumAccessOpenMS input , list output , list extraction_coordinates , double mz_extraction_window ,  ppm , bytes filter ):
        assert isinstance(input, SpectrumAccessOpenMS), 'arg input wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, OSChromatogram) for elemt_rec in output), 'arg output wrong type'
        assert isinstance(extraction_coordinates, list) and all(isinstance(elemt_rec, ExtractionCoordinates) for elemt_rec in extraction_coordinates), 'arg extraction_coordinates wrong type'
        assert isinstance(mz_extraction_window, float), 'arg mz_extraction_window wrong type'
        assert isinstance(ppm, (int, long)), 'arg ppm wrong type'
        assert isinstance(filter, bytes), 'arg filter wrong type'

        # cdef shared_ptr[_SpectrumAccessOpenMS] input_input = input.inst
        cdef libcpp_vector[shared_ptr[_OSChromatogram]] v1
        cdef OSChromatogram output_rec
        for output_rec in output:
            v1.push_back(output_rec.inst)
        cdef libcpp_vector[_ExtractionCoordinates] * v2 = new libcpp_vector[_ExtractionCoordinates]()
        cdef ExtractionCoordinates item2
        for item2 in extraction_coordinates:
            v2.push_back(deref(item2.inst.get()))
    
        # Call 
        self.inst.get().extractChromatograms((<shared_ptr[_ISpectrumAccess]>input.inst), v1, deref(v2), (<double>mz_extraction_window), (<bool>ppm), (_String(<char *>filter)))

        del v2

        # gather results
        replace = list()
        cdef libcpp_vector[shared_ptr[_OSChromatogram]].iterator it_output = v1.begin()
        while it_output != v1.end():
           output_rec = OSChromatogram.__new__(OSChromatogram)
           output_rec.inst = deref(it_output)
           replace.append(output_rec)
           inc(it_output)
        # replace old vector with new contents
        output[:] = []
        for output_rec_b in replace:
            output.append(output_rec_b)

    def _extractChromatograms_OpenMSCached(self, SpectrumAccessOpenMSCached input , list output , list extraction_coordinates , double mz_extraction_window ,  ppm , bytes filter ):
        assert isinstance(input, SpectrumAccessOpenMSCached), 'arg input wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, OSChromatogram) for elemt_rec in output), 'arg output wrong type'
        assert isinstance(extraction_coordinates, list) and all(isinstance(elemt_rec, ExtractionCoordinates) for elemt_rec in extraction_coordinates), 'arg extraction_coordinates wrong type'
        assert isinstance(mz_extraction_window, float), 'arg mz_extraction_window wrong type'
        assert isinstance(ppm, (int, long)), 'arg ppm wrong type'
        assert isinstance(filter, bytes), 'arg filter wrong type'

        # cdef shared_ptr[_SpectrumAccessOpenMSCached] input_input = input.inst
        cdef libcpp_vector[shared_ptr[_OSChromatogram]] v1
        cdef OSChromatogram output_rec
        for output_rec in output:
            v1.push_back(output_rec.inst)
        cdef libcpp_vector[_ExtractionCoordinates] * v2 = new libcpp_vector[_ExtractionCoordinates]()
        cdef ExtractionCoordinates item2
        for item2 in extraction_coordinates:
            v2.push_back(deref(item2.inst.get()))
    
        # Call 
        self.inst.get().extractChromatograms((<shared_ptr[_ISpectrumAccess]>input.inst), v1, deref(v2), (<double>mz_extraction_window), (<bool>ppm), (_String(<char *>filter)))

        del v2

        # gather results
        replace = list()
        cdef libcpp_vector[shared_ptr[_OSChromatogram]].iterator it_output = v1.begin()
        while it_output != v1.end():
           output_rec = OSChromatogram.__new__(OSChromatogram)
           output_rec.inst = deref(it_output)
           replace.append(output_rec)
           inc(it_output)
        # replace old vector with new contents
        output[:] = []
        for output_rec_b in replace:
            output.append(output_rec_b)


