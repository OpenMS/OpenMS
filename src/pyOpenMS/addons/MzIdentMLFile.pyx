


    def load(self, filename, protein_ids, peptide_ids):
        """
        load(self, filename, protein_ids, peptide_ids)

        Loads the identifications from a MzIdentML file.

        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for the peptide_ids parameter.

        :param filename: Path to the MzIdentML file to load
        :type filename: str or bytes or String
        :param protein_ids: List to store protein identifications (modified in place)
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (deprecated) and
                           PeptideIdentificationList (recommended).
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList

        Example::

            # Recommended: Use PeptideIdentificationList
            protein_ids = []
            peptide_ids = pyopenms.PeptideIdentificationList()
            pyopenms.MzIdentMLFile().load("test.mzid", protein_ids, peptide_ids)

            # Backward compatible: Python list (deprecated)
            protein_ids = []
            peptide_ids = []
            pyopenms.MzIdentMLFile().load("test.mzid", protein_ids, peptide_ids)

        .. deprecated:: 3.5
            Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5.
            Use PeptideIdentificationList instead. The list interface will be removed
            in a future version.
        """
        import warnings

        # Convert protein_ids list to C++ vector
        cdef libcpp_vector[_ProteinIdentification] * c_protein_ids = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification prot_item
        for prot_item in protein_ids:
            c_protein_ids.push_back(deref(prot_item.inst.get()))

        # Check if peptide_ids is a Python list (deprecated) or PeptideIdentificationList (new API)
        cdef PeptideIdentificationList temp_peptide_ids
        cdef PeptideIdentificationList pep_id_list

        if isinstance(peptide_ids, list):
            # Deprecated: Python list interface
            warnings.warn(
                "Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5. "
                "Use PeptideIdentificationList instead: peptide_ids = pyopenms.PeptideIdentificationList()",
                DeprecationWarning,
                stacklevel=2
            )
            temp_peptide_ids = PeptideIdentificationList()
            temp_peptide_ids.extend(peptide_ids)

            # Call the C++ load method
            self.inst.get().load(
                deref(convString(filename).get()),
                deref(c_protein_ids),
                deref(temp_peptide_ids.inst.get())
            )

            # Copy results back to Python list
            peptide_ids[:] = list(temp_peptide_ids)
        else:
            # New API: Use PeptideIdentificationList directly
            pep_id_list = <PeptideIdentificationList>peptide_ids
            self.inst.get().load(
                deref(convString(filename).get()),
                deref(c_protein_ids),
                deref(pep_id_list.inst.get())
            )

        # Copy protein_ids back to Python list
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot = c_protein_ids.begin()
        cdef list replace_prot = []
        while it_prot != c_protein_ids.end():
            prot_item = ProteinIdentification.__new__(ProteinIdentification)
            prot_item.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot)))
            replace_prot.append(prot_item)
            inc(it_prot)
        protein_ids[:] = replace_prot
        del c_protein_ids

    def store(self, filename, protein_ids, peptide_ids):
        """
        store(self, filename, protein_ids, peptide_ids)

        Stores the identifications to a MzIdentML file.

        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for the peptide_ids parameter.

        :param filename: Path to the MzIdentML file to store
        :type filename: str or bytes or String
        :param protein_ids: List of protein identifications to store
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container of peptide identifications to store.
                           Accepts both list[PeptideIdentification] (deprecated) and
                           PeptideIdentificationList (recommended).
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList

        Example::

            # Recommended: Use PeptideIdentificationList
            protein_ids = [...]
            peptide_ids = pyopenms.PeptideIdentificationList()
            peptide_ids.extend([...])
            pyopenms.MzIdentMLFile().store("test.mzid", protein_ids, peptide_ids)

            # Backward compatible: Python list (deprecated)
            protein_ids = [...]
            peptide_ids = [...]
            pyopenms.MzIdentMLFile().store("test.mzid", protein_ids, peptide_ids)

        .. deprecated:: 3.5
            Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5.
            Use PeptideIdentificationList instead. The list interface will be removed
            in a future version.
        """
        import warnings

        # Convert protein_ids list to C++ vector
        cdef libcpp_vector[_ProteinIdentification] * c_protein_ids = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification prot_item
        for prot_item in protein_ids:
            c_protein_ids.push_back(deref(prot_item.inst.get()))

        # Check if peptide_ids is a Python list (deprecated) or PeptideIdentificationList (new API)
        cdef PeptideIdentificationList temp_peptide_ids
        cdef PeptideIdentificationList pep_id_list

        if isinstance(peptide_ids, list):
            # Deprecated: Python list interface
            warnings.warn(
                "Passing a Python list for peptide_ids is deprecated since pyOpenMS 3.5. "
                "Use PeptideIdentificationList instead: peptide_ids = pyopenms.PeptideIdentificationList()",
                DeprecationWarning,
                stacklevel=2
            )
            temp_peptide_ids = PeptideIdentificationList()
            temp_peptide_ids.extend(peptide_ids)

            self.inst.get().store(
                deref(convString(filename).get()),
                deref(c_protein_ids),
                deref(temp_peptide_ids.inst.get())
            )
        else:
            # New API: Use PeptideIdentificationList directly
            pep_id_list = <PeptideIdentificationList>peptide_ids
            self.inst.get().store(
                deref(convString(filename).get()),
                deref(c_protein_ids),
                deref(pep_id_list.inst.get())
            )

        del c_protein_ids

