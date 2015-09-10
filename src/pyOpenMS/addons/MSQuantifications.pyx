


    def registerExperiment(self, MSExperiment exp, list labels):

        assert isinstance(labels, list) and all(isinstance(li, list) and all( 
                isinstance(inner_li, list) and all(
                    isinstance(inner_pair, list) and len(inner_pair) == 2 and
                        isinstance(inner_pair[0], String) and isinstance(inner_pair[1], float)
                for inner_pair in inner_li)
            for inner_li in li) 
        for li in labels), 'arg proteins wrong type'

        cdef libcpp_vector[libcpp_vector[libcpp_pair[_String, double]]] * v0 = new libcpp_vector[libcpp_vector[libcpp_pair[_String, double]]]()
        cdef libcpp_vector[libcpp_pair[_String, double] ] * v0_rec = new libcpp_vector[libcpp_pair[_String, double]]()
        cdef String innerstring_
        cdef libcpp_vector[libcpp_pair[_String, double] ].iterator it_labels_rec
        cdef libcpp_pair[_String, double] * aPair
        for labels_rec in labels:
            v0_rec.clear()
            for item0_rec in labels_rec:
                
                # Build new pair, push it, delete it 
                innerstring_ = item0_rec.first
                aPair = new libcpp_pair[ _String, double]( deref(innerstring_.inst.get()), <double>item0_rec.second )
                v0_rec.push_back(deref(aPair))
                del aPair

            v0.push_back(deref(v0_rec))

        self.inst.get().registerExperiment( (deref( exp.inst.get() ))  , (deref(v0)) )

        del v0
        del v0_rec

    # def getRatios(self):
    #     cdef libcpp_map[_String, _Ratio] res 
    #     c_res = self.inst.get().getRatios()
    #     #
    #     ## Get the data back from C++
    #     #
    #     replace = dict()
    #     cdef libcpp_map[_String, _Ratio].iterator it_ripped = c_res.begin()
    #     cdef _Ratio myRatio
    #     cdef Ratio py_ratio
    #     while it_ripped != c_res.end():
    #         
    #         myRatio = deref(it_ripped).second
    #         py_ratio = Ratio.__new__(Ratio)
    #         py_ratio.inst = shared_ptr[_Ratio](new _Ratio(myRatio))

    #         replace[ <libcpp_string>deref(it_ripped).first ]  = py_ratio
    #         inc(it_ripped)

    #     return replace



