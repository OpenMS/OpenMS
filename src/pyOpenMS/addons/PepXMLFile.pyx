


    def load(self, filename, protein_ids, peptide_ids, experiment_name=None, lookup=None):
        """
        load(self, filename, protein_ids, peptide_ids, experiment_name=None, lookup=None)

        Loads the identifications from a PepXML file.

        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for the peptide_ids parameter.

        :param filename: Path to the PepXML file to load
        :type filename: str or bytes or String
        :param protein_ids: List to store protein identifications (modified in place)
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (deprecated) and
                           PeptideIdentificationList (recommended).
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param experiment_name: Optional experiment name
        :type experiment_name: str or bytes or String, optional
        :param lookup: Optional spectrum metadata lookup structure
        :type lookup: SpectrumMetaDataLookup, optional

        Example::

            # Recommended: Use PeptideIdentificationList
            protein_ids = []
            peptide_ids = pyopenms.PeptideIdentificationList()
            pyopenms.PepXMLFile().load("test.pep.xml", protein_ids, peptide_ids)

            # Backward compatible: Python list (deprecated)
            protein_ids = []
            peptide_ids = []
            pyopenms.PepXMLFile().load("test.pep.xml", protein_ids, peptide_ids)

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
        cdef SpectrumMetaDataLookup c_lookup

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

            # Call the appropriate C++ method based on arguments
            if experiment_name is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(temp_peptide_ids.inst.get())
                )
            elif lookup is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(temp_peptide_ids.inst.get()),
                    deref(convString(experiment_name).get())
                )
            else:
                c_lookup = <SpectrumMetaDataLookup>lookup
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(temp_peptide_ids.inst.get()),
                    deref(convString(experiment_name).get()),
                    deref(c_lookup.inst.get())
                )

            # Copy results back to Python list
            peptide_ids[:] = list(temp_peptide_ids)
        else:
            # New API: Use PeptideIdentificationList directly
            pep_id_list = <PeptideIdentificationList>peptide_ids
            if experiment_name is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(pep_id_list.inst.get())
                )
            elif lookup is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(pep_id_list.inst.get()),
                    deref(convString(experiment_name).get())
                )
            else:
                c_lookup = <SpectrumMetaDataLookup>lookup
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(pep_id_list.inst.get()),
                    deref(convString(experiment_name).get()),
                    deref(c_lookup.inst.get())
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

    def store(self, filename, protein_ids, peptide_ids, mz_file=None, mz_name=None,
              peptideprophet_analyzed=False, rt_tolerance=0.0):
        """
        store(self, filename, protein_ids, peptide_ids, mz_file=None, mz_name=None,
              peptideprophet_analyzed=False, rt_tolerance=0.0)

        Stores the identifications to a PepXML file.

        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for the peptide_ids parameter.

        :param filename: Path to the PepXML file to store
        :type filename: str or bytes or String
        :param protein_ids: List of protein identifications to store
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container of peptide identifications to store.
                           Accepts both list[PeptideIdentification] (deprecated) and
                           PeptideIdentificationList (recommended).
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param mz_file: Optional mz file name
        :type mz_file: str or bytes or String, optional
        :param mz_name: Optional mz name
        :type mz_name: str or bytes or String, optional
        :param peptideprophet_analyzed: Whether data was analyzed with PeptideProphet
        :type peptideprophet_analyzed: bool, optional
        :param rt_tolerance: RT tolerance
        :type rt_tolerance: float, optional

        Example::

            # Recommended: Use PeptideIdentificationList
            protein_ids = [...]
            peptide_ids = pyopenms.PeptideIdentificationList()
            peptide_ids.extend([...])
            pyopenms.PepXMLFile().store("test.pep.xml", protein_ids, peptide_ids)

            # Backward compatible: Python list (deprecated)
            protein_ids = [...]
            peptide_ids = [...]
            pyopenms.PepXMLFile().store("test.pep.xml", protein_ids, peptide_ids)

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

            if mz_file is None:
                self.inst.get().store(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(temp_peptide_ids.inst.get())
                )
            else:
                self.inst.get().store(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(temp_peptide_ids.inst.get()),
                    deref(convString(mz_file).get()),
                    deref(convString(mz_name).get()),
                    <bint>peptideprophet_analyzed,
                    <double>rt_tolerance
                )
        else:
            # New API: Use PeptideIdentificationList directly
            pep_id_list = <PeptideIdentificationList>peptide_ids
            if mz_file is None:
                self.inst.get().store(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(pep_id_list.inst.get())
                )
            else:
                self.inst.get().store(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(pep_id_list.inst.get()),
                    deref(convString(mz_file).get()),
                    deref(convString(mz_name).get()),
                    <bint>peptideprophet_analyzed,
                    <double>rt_tolerance
                )

        del c_protein_ids

