

    def rip(self, ripped, list proteins, list peptides):

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

        self.inst.get().rip(c_ripped, deref(c_protein_vec), deref(c_peptide_vec))
        del c_peptide_vec
        del c_protein_vec
        del c_protein_vec_inner
        del c_peptide_vec_inner

