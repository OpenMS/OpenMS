import unittest
import os
import tempfile
import warnings
import pyopenms

class TestBackwardCompatibilityIdXMLFile(unittest.TestCase):
    """Test backward compatibility for IdXMLFile with Python lists"""

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.idXML").encode()

    def test_load_with_python_list(self):
        """Test that IdXMLFile.load() works with Python lists (deprecated API)"""
        idxml_file = pyopenms.IdXMLFile()
        protein_ids = []
        peptide_ids = []  # Deprecated: Python list

        # This should work but emit a DeprecationWarning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            idxml_file.load(self.filename, protein_ids, peptide_ids)

            # Verify deprecation warning was emitted
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[0].category, DeprecationWarning))
            self.assertIn("deprecated", str(w[0].message).lower())

        # Verify the results
        self.assertEqual(len(protein_ids), 1)
        self.assertEqual(len(peptide_ids), 3)
        self.assertIsInstance(peptide_ids, list)
        self.assertIsInstance(peptide_ids[0], pyopenms.PeptideIdentification)

    def test_load_with_peptide_identification_list(self):
        """Test that IdXMLFile.load() works with PeptideIdentificationList (new API)"""
        idxml_file = pyopenms.IdXMLFile()
        protein_ids = []
        peptide_ids = pyopenms.PeptideIdentificationList()  # New API

        # This should work without any warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            idxml_file.load(self.filename, protein_ids, peptide_ids)

            # Verify no deprecation warning was emitted
            deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            self.assertEqual(len(deprecation_warnings), 0)

        # Verify the results
        self.assertEqual(len(protein_ids), 1)
        self.assertEqual(peptide_ids.size(), 3)
        self.assertIsInstance(peptide_ids, pyopenms.PeptideIdentificationList)

    def test_store_with_python_list(self):
        """Test that IdXMLFile.store() works with Python lists (deprecated API)"""
        # First load data using new API
        idxml_file = pyopenms.IdXMLFile()
        protein_ids = []
        peptide_ids = pyopenms.PeptideIdentificationList()
        idxml_file.load(self.filename, protein_ids, peptide_ids)

        # Convert to list for deprecated API test
        peptide_ids_list = list(peptide_ids)

        # Now store using deprecated Python list
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.idXML') as tmpfile:
            temp_filename = tmpfile.name

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                idxml_file.store(temp_filename.encode(), protein_ids, peptide_ids_list)

                # Verify deprecation warning was emitted
                self.assertEqual(len(w), 1)
                self.assertTrue(issubclass(w[0].category, DeprecationWarning))

            # Verify we can read it back
            protein_ids_2 = []
            peptide_ids_2 = pyopenms.PeptideIdentificationList()
            idxml_file.load(temp_filename.encode(), protein_ids_2, peptide_ids_2)

            self.assertEqual(len(protein_ids_2), 1)
            self.assertEqual(peptide_ids_2.size(), 3)
        finally:
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)

    def test_store_with_peptide_identification_list(self):
        """Test that IdXMLFile.store() works with PeptideIdentificationList (new API)"""
        # First load data
        idxml_file = pyopenms.IdXMLFile()
        protein_ids = []
        peptide_ids = pyopenms.PeptideIdentificationList()
        idxml_file.load(self.filename, protein_ids, peptide_ids)

        # Now store using new-style PeptideIdentificationList
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.idXML') as tmpfile:
            temp_filename = tmpfile.name

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                idxml_file.store(temp_filename.encode(), protein_ids, peptide_ids)

                # Verify no deprecation warning was emitted
                deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
                self.assertEqual(len(deprecation_warnings), 0)

            # Verify we can read it back
            protein_ids_2 = []
            peptide_ids_2 = pyopenms.PeptideIdentificationList()
            idxml_file.load(temp_filename.encode(), protein_ids_2, peptide_ids_2)

            self.assertEqual(len(protein_ids_2), 1)
            self.assertEqual(peptide_ids_2.size(), 3)
        finally:
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)


class TestBackwardCompatibilityPepXMLFile(unittest.TestCase):
    """Test backward compatibility for PepXMLFile with Python lists"""

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.pep.xml").encode()

    def test_load_with_python_list(self):
        """Test that PepXMLFile.load() works with Python lists (deprecated API)"""
        pepxml_file = pyopenms.PepXMLFile()
        protein_ids = []
        peptide_ids = []  # Deprecated: Python list

        # This should work but emit a DeprecationWarning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            pepxml_file.load(self.filename, protein_ids, peptide_ids)

            # Verify deprecation warning was emitted
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[0].category, DeprecationWarning))

        # Verify the results
        self.assertEqual(len(protein_ids), 3)
        self.assertEqual(len(peptide_ids), 19)
        self.assertIsInstance(peptide_ids, list)
        self.assertIsInstance(peptide_ids[0], pyopenms.PeptideIdentification)

    def test_load_with_peptide_identification_list(self):
        """Test that PepXMLFile.load() works with PeptideIdentificationList (new API)"""
        pepxml_file = pyopenms.PepXMLFile()
        protein_ids = []
        peptide_ids = pyopenms.PeptideIdentificationList()  # New API

        # This should work without any warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            pepxml_file.load(self.filename, protein_ids, peptide_ids)

            # Verify no deprecation warning was emitted
            deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            self.assertEqual(len(deprecation_warnings), 0)

        # Verify the results
        self.assertEqual(len(protein_ids), 3)
        self.assertEqual(peptide_ids.size(), 19)
        self.assertIsInstance(peptide_ids, pyopenms.PeptideIdentificationList)

    def test_store_with_python_list(self):
        """Test that PepXMLFile.store() works with Python lists (deprecated API)"""
        # First load data using new API
        pepxml_file = pyopenms.PepXMLFile()
        protein_ids = []
        peptide_ids = pyopenms.PeptideIdentificationList()
        pepxml_file.load(self.filename, protein_ids, peptide_ids)

        # Convert to list for deprecated API test
        peptide_ids_list = list(peptide_ids)

        # Now store using deprecated Python list
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.pep.xml') as tmpfile:
            temp_filename = tmpfile.name

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                pepxml_file.store(temp_filename.encode(), protein_ids, peptide_ids_list)

                # Verify deprecation warning was emitted
                self.assertEqual(len(w), 1)
                self.assertTrue(issubclass(w[0].category, DeprecationWarning))

            # Verify we can read it back
            protein_ids_2 = []
            peptide_ids_2 = pyopenms.PeptideIdentificationList()
            pepxml_file.load(temp_filename.encode(), protein_ids_2, peptide_ids_2)

            # Note: PepXML format has limitations and may not preserve all protein IDs
            # during roundtrip, but peptide IDs should be preserved
            self.assertGreaterEqual(len(protein_ids_2), 1)
            self.assertEqual(peptide_ids_2.size(), 19)
        finally:
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)


class TestBackwardCompatibilityMzIdentMLFile(unittest.TestCase):
    """Test backward compatibility for MzIdentMLFile with Python lists"""

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.mzid").encode()

    def test_load_with_python_list(self):
        """Test that MzIdentMLFile.load() works with Python lists (deprecated API)"""
        mzid_file = pyopenms.MzIdentMLFile()
        protein_ids = []
        peptide_ids = []  # Deprecated: Python list

        # This should work but emit a DeprecationWarning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mzid_file.load(self.filename, protein_ids, peptide_ids)

            # Verify deprecation warning was emitted
            self.assertEqual(len(w), 1)
            self.assertTrue(issubclass(w[0].category, DeprecationWarning))

        # Verify the results
        self.assertGreater(len(protein_ids), 0)
        self.assertGreater(len(peptide_ids), 0)
        self.assertIsInstance(peptide_ids, list)
        self.assertIsInstance(peptide_ids[0], pyopenms.PeptideIdentification)

    def test_load_with_peptide_identification_list(self):
        """Test that MzIdentMLFile.load() works with PeptideIdentificationList (new API)"""
        mzid_file = pyopenms.MzIdentMLFile()
        protein_ids = []
        peptide_ids = pyopenms.PeptideIdentificationList()  # New API

        # This should work without any warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mzid_file.load(self.filename, protein_ids, peptide_ids)

            # Verify no deprecation warning was emitted
            deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            self.assertEqual(len(deprecation_warnings), 0)

        # Verify the results
        self.assertGreater(len(protein_ids), 0)
        self.assertGreater(peptide_ids.size(), 0)
        self.assertIsInstance(peptide_ids, pyopenms.PeptideIdentificationList)

    def test_store_with_python_list(self):
        """Test that MzIdentMLFile.store() works with Python lists (deprecated API)"""
        # First load data using new API
        mzid_file = pyopenms.MzIdentMLFile()
        protein_ids = []
        peptide_ids = pyopenms.PeptideIdentificationList()
        mzid_file.load(self.filename, protein_ids, peptide_ids)

        original_peptide_count = peptide_ids.size()

        # Convert to list for deprecated API test
        peptide_ids_list = list(peptide_ids)

        # Now store using deprecated Python list
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.mzid') as tmpfile:
            temp_filename = tmpfile.name

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                mzid_file.store(temp_filename.encode(), protein_ids, peptide_ids_list)

                # Verify deprecation warning was emitted
                self.assertEqual(len(w), 1)
                self.assertTrue(issubclass(w[0].category, DeprecationWarning))

            # Verify we can read it back
            protein_ids_2 = []
            peptide_ids_2 = pyopenms.PeptideIdentificationList()
            mzid_file.load(temp_filename.encode(), protein_ids_2, peptide_ids_2)

            self.assertGreater(len(protein_ids_2), 0)
            self.assertEqual(peptide_ids_2.size(), original_peptide_count)
        finally:
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)

    def test_store_with_peptide_identification_list(self):
        """Test that MzIdentMLFile.store() works with PeptideIdentificationList (new API)"""
        # First load data
        mzid_file = pyopenms.MzIdentMLFile()
        protein_ids = []
        peptide_ids = pyopenms.PeptideIdentificationList()
        mzid_file.load(self.filename, protein_ids, peptide_ids)

        original_peptide_count = peptide_ids.size()

        # Now store using new-style PeptideIdentificationList
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.mzid') as tmpfile:
            temp_filename = tmpfile.name

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                mzid_file.store(temp_filename.encode(), protein_ids, peptide_ids)

                # Verify no deprecation warning was emitted
                deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
                self.assertEqual(len(deprecation_warnings), 0)

            # Verify we can read it back
            protein_ids_2 = []
            peptide_ids_2 = pyopenms.PeptideIdentificationList()
            mzid_file.load(temp_filename.encode(), protein_ids_2, peptide_ids_2)

            self.assertGreater(len(protein_ids_2), 0)
            self.assertEqual(peptide_ids_2.size(), original_peptide_count)
        finally:
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)


if __name__ == '__main__':
    unittest.main()
