

    def rip(self, ripped, list proteins, list peptides):
 
        # TODO assert correct type of ripped
        assert isinstance(proteins, list) and all(isinstance(li, ProteinIdentification) for li in proteins), 'arg proteins wrong type'
        cdef libcpp_vector[_ProteinIdentification] * c_protein_vec = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in proteins:
           c_protein_vec.push_back(deref(item1.inst.get()))
 
        assert isinstance(peptides, list) and all(isinstance(li, PeptideIdentification) for li in peptides), 'arg peptides wrong type'
        cdef libcpp_vector[_PeptideIdentification] * c_peptide_vec = new libcpp_vector[_PeptideIdentification]()
        cdef PeptideIdentification item0
        for item0 in peptides:
           c_peptide_vec.push_back(deref(item0.inst.get()))
 
        cdef libcpp_map[_String, libcpp_pair[ libcpp_vector[_ProteinIdentification], libcpp_vector[_PeptideIdentification]]] c_ripped
        # declaration for loop 
        cdef libcpp_vector[_ProteinIdentification] * c_protein_vec_inner = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification i_item1
        cdef libcpp_vector[_PeptideIdentification] * c_peptide_vec_inner = new libcpp_vector[_PeptideIdentification]()
        cdef PeptideIdentification i_item0
        cdef libcpp_pair[ libcpp_vector[_ProteinIdentification], libcpp_vector[_PeptideIdentification]] * aPair
        for k,v in ripped.iteritems():
 
            c_protein_vec_inner.clear()
            c_peptide_vec_inner.clear()
 
            prot_vec = v[0]
            assert isinstance(prot_vec, list) and all(isinstance(li, ProteinIdentification) for li in prot_vec), 'arg proteins wrong type'
            for i_item1 in prot_vec:
               c_protein_vec_inner.push_back(deref(i_item1.inst.get()))
 
            pep_vec = v[1]
            assert isinstance(pep_vec, list) and all(isinstance(li, PeptideIdentification) for li in pep_vec), 'arg peptides wrong type'
            for i_item0 in pep_vec:
               c_peptide_vec_inner.push_back(deref(i_item0.inst.get()))
 
            aPair = new libcpp_pair[ libcpp_vector[_ProteinIdentification], libcpp_vector[_PeptideIdentification]]( deref(c_protein_vec_inner), deref(c_peptide_vec_inner) )
 
            assert isinstance(k, bytes), 'arg key in label_identifiers wrong type'
            assert isinstance(v, float), 'arg value in label_identifiers wrong type'
            c_ripped[ _String(<char *>k) ] = deref(aPair)
 
            del aPair
 
        # call
        self.inst.get().rip(c_ripped, deref(c_protein_vec), deref(c_peptide_vec))

        # get the data back from C++
        replace = dict()
        cdef libcpp_map[_String, libcpp_pair[ libcpp_vector[_ProteinIdentification], libcpp_vector[_PeptideIdentification]]].iterator it_ripped = c_ripped.begin()
        cdef libcpp_pair[ libcpp_vector[_ProteinIdentification], libcpp_vector[_PeptideIdentification]] anotherPair
        cdef libcpp_vector[_ProteinIdentification] another_c_protein_vec_inner 
        cdef libcpp_vector[_PeptideIdentification] another_c_peptide_vec_inner 
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot
        cdef libcpp_vector[_PeptideIdentification].iterator it_pep
        cdef ProteinIdentification item_py_result_prot
        cdef PeptideIdentification item_py_result_pep
        while it_ripped != c_ripped.end():
            # get the pair and the two vectors
            anotherPair = deref(it_ripped).second
            another_c_protein_vec_inner = anotherPair.first
            another_c_peptide_vec_inner = anotherPair.second

            replace_inner_protein = []
            it_prot = another_c_protein_vec_inner.begin()
            while it_prot != another_c_protein_vec_inner.end():
                item_py_result_prot = ProteinIdentification.__new__(ProteinIdentification)   
                item_py_result_prot.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot)))
                replace_inner_protein.append(item_py_result_prot)
                inc(it_prot)

            replace_inner_peptide = []
            it_pep = another_c_peptide_vec_inner.begin()
            while it_pep != another_c_peptide_vec_inner.end():
                item_py_result_pep = PeptideIdentification.__new__(PeptideIdentification)   
                item_py_result_pep.inst = shared_ptr[_PeptideIdentification](new _PeptideIdentification(deref(it_pep)))
                replace_inner_peptide.append(item_py_result_pep)
                inc(it_pep)

            replace[ <libcpp_string>deref(it_ripped).first ]  = [ replace_inner_protein, replace_inner_peptide] 
            inc(it_ripped)
        ripped.clear()
        ripped.update(replace)

        del c_peptide_vec
        del c_protein_vec
        del c_protein_vec_inner
        del c_peptide_vec_inner
 

