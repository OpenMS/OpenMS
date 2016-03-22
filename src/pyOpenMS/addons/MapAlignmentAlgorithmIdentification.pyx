


    def align(self, list ids , list trafos, int ref_index):
        assert isinstance(ids, list) and all(isinstance(li, list) and all(isinstance(li_, PeptideIdentification) for li_ in li) for li in ids), 'arg ids wrong type'
        assert isinstance(trafos, list) and all(isinstance(li, TransformationDescription) for li in trafos), 'arg trafos wrong type'
        cdef libcpp_vector[libcpp_vector[_PeptideIdentification]] * v0 = new libcpp_vector[libcpp_vector[_PeptideIdentification]]()
        cdef libcpp_vector[_PeptideIdentification] * v0_rec = new libcpp_vector[_PeptideIdentification]()
        cdef PeptideIdentification item0_rec
        cdef libcpp_vector[_PeptideIdentification].iterator it_ids_rec
        for ids_rec in ids:
            v0_rec.clear()
            for item0_rec in ids_rec:
                v0_rec.push_back(deref(item0_rec.inst.get()))
            v0.push_back(deref(v0_rec))
        cdef libcpp_vector[_TransformationDescription] * v1 = new libcpp_vector[_TransformationDescription]()
        cdef TransformationDescription item1
        for item1 in trafos:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().align(deref(v0), deref(v1), ref_index)
        cdef libcpp_vector[_TransformationDescription].iterator it_trafos = v1.begin()
        replace_0 = []
        while it_trafos != v1.end():
            item1 = TransformationDescription.__new__(TransformationDescription)
            item1.inst = shared_ptr[_TransformationDescription](new _TransformationDescription(deref(it_trafos)))
            replace_0.append(item1)
            inc(it_trafos)
        trafos[:] = replace_0
        del v1
        cdef libcpp_vector[libcpp_vector[_PeptideIdentification]].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            it_ids_rec = deref(it_ids).begin()
            replace_1 = []
            while it_ids_rec != deref(it_ids).end():
                item0_rec = PeptideIdentification.__new__(PeptideIdentification)
                item0_rec.inst = shared_ptr[_PeptideIdentification](new _PeptideIdentification(deref(it_ids_rec)))
                replace_1.append(item0_rec)
                inc(it_ids_rec)
            replace_0.append(replace_1)
            inc(it_ids)
        ids[:] = replace_0
        del v0

