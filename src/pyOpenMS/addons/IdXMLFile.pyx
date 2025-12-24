


    def load(self, filename, protein_ids, peptide_ids, document_id=None):
        """
        load(self, filename, protein_ids, peptide_ids, document_id=None)

        Loads the identifications from an idXML file.

        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.

        :param filename: Path to the idXML file to load
        :type filename: str or bytes or String
        :param protein_ids: List to store protein identifications (modified in place)
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container to store peptide identifications (modified in place).
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param document_id: Optional document identifier (output parameter)
        :type document_id: String, optional

        Example::

            # New style (3.5+)
            protein_ids = []
            peptide_ids = pyopenms.PeptideIdentificationList()
            pyopenms.IdXMLFile().load("test.idXML", protein_ids, peptide_ids)

            # Old style (backward compatible)
            protein_ids = []
            peptide_ids = []
            pyopenms.IdXMLFile().load("test.idXML", protein_ids, peptide_ids)
        """
        # Convert protein_ids list to C++ vector
        cdef libcpp_vector[_ProteinIdentification] * c_protein_ids = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification prot_item
        for prot_item in protein_ids:
            c_protein_ids.push_back(deref(prot_item.inst.get()))

        # Check if peptide_ids is a Python list (old API) or PeptideIdentificationList (new API)
        cdef PeptideIdentificationList temp_peptide_ids
        cdef PeptideIdentificationList pep_id_list

        if isinstance(peptide_ids, list):
            # Old API: Convert Python list to PeptideIdentificationList
            temp_peptide_ids = PeptideIdentificationList()
            temp_peptide_ids.extend(peptide_ids)

            # Call the C++ load method
            if document_id is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(temp_peptide_ids.inst.get())
                )
            else:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(temp_peptide_ids.inst.get()),
                    deref(convString(document_id).get())
                )

            # Copy results back to Python list
            peptide_ids[:] = list(temp_peptide_ids)
        else:
            # New API: Use PeptideIdentificationList directly
            pep_id_list = <PeptideIdentificationList>peptide_ids
            if document_id is None:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(pep_id_list.inst.get())
                )
            else:
                self.inst.get().load(
                    deref(convString(filename).get()),
                    deref(c_protein_ids),
                    deref(pep_id_list.inst.get()),
                    deref(convString(document_id).get())
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

    def store(self, filename, protein_ids, peptide_ids, document_id=b""):
        """
        store(self, filename, protein_ids, peptide_ids, document_id="")

        Stores the identifications to an idXML file.

        Provides backward compatibility by accepting both Python lists and
        PeptideIdentificationList objects for peptide_ids parameter.

        :param filename: Path to the idXML file to store
        :type filename: str or bytes or String
        :param protein_ids: List of protein identifications to store
        :type protein_ids: list[ProteinIdentification]
        :param peptide_ids: Container of peptide identifications to store.
                           Accepts both list[PeptideIdentification] (pre-3.5 style) and
                           PeptideIdentificationList (3.5+ style)
        :type peptide_ids: list[PeptideIdentification] or PeptideIdentificationList
        :param document_id: Optional document identifier
        :type document_id: str or bytes or String, optional

        Example::

            # New style (3.5+)
            protein_ids = [...]
            peptide_ids = pyopenms.PeptideIdentificationList()
            peptide_ids.extend([...])
            pyopenms.IdXMLFile().store("test.idXML", protein_ids, peptide_ids)

            # Old style (backward compatible)
            protein_ids = [...]
            peptide_ids = [...]
            pyopenms.IdXMLFile().store("test.idXML", protein_ids, peptide_ids)
        """
        # Convert protein_ids list to C++ vector
        cdef libcpp_vector[_ProteinIdentification] * c_protein_ids = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification prot_item
        for prot_item in protein_ids:
            c_protein_ids.push_back(deref(prot_item.inst.get()))

        # Check if peptide_ids is a Python list (old API) or PeptideIdentificationList (new API)
        cdef PeptideIdentificationList temp_peptide_ids
        cdef PeptideIdentificationList pep_id_list

        if isinstance(peptide_ids, list):
            # Old API: Convert Python list to PeptideIdentificationList
            temp_peptide_ids = PeptideIdentificationList()
            temp_peptide_ids.extend(peptide_ids)

            self.inst.get().store(
                deref(convString(filename).get()),
                deref(c_protein_ids),
                deref(temp_peptide_ids.inst.get()),
                deref(convString(document_id).get())
            )
        else:
            # New API: Use PeptideIdentificationList directly
            pep_id_list = <PeptideIdentificationList>peptide_ids
            self.inst.get().store(
                deref(convString(filename).get()),
                deref(c_protein_ids),
                deref(pep_id_list.inst.get()),
                deref(convString(document_id).get())
            )

        del c_protein_ids

