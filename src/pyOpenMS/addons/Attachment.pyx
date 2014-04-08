


    property tableRows:
    
        def __set__(self, list tableRows):
            # libcpp_vector[ libcpp_vector[String] ]
            cdef libcpp_vector[libcpp_vector[_String]] * v0 = new libcpp_vector[libcpp_vector[_String]]()
            cdef libcpp_vector[_String] * v0_rec = new libcpp_vector[_String]()
            cdef bytes item0
            for tableRows_rec in tableRows:
                v0_rec.clear()
                for item0 in tableRows_rec:
                    v0_rec.push_back(_String(<char *>item0))
                v0.push_back(deref(v0_rec))
            self.inst.get().tableRows = deref(v0)
            del v0
            del v0_rec
    

        def __get__(self):
            c_res = self.inst.get().tableRows
            # libcpp_vector[ libcpp_vector[String] ]
            cdef replace = []
            cdef libcpp_vector[libcpp_vector[_String] ].iterator it = c_res.begin()
            cdef libcpp_vector[_String ].iterator it_inner
            while it != c_res.end():
                replace_inner = []
                it_inner = deref(it).begin()
                while it_inner != deref(it).end():
                    replace_inner.append(<char*>deref(it_inner).c_str())
                    inc(it_inner)
                replace.append(replace_inner)
                inc(it)
            return replace 
