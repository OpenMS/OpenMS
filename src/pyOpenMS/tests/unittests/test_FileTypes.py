"""
Test FileTypes class wrapping.
Tests that:
1. Static methods are properly marked with @staticmethod
2. All FileType enum values are accessible
3. FileProperties enum class is properly wrapped
4. FilterLayout enum class is properly wrapped
5. Static methods work correctly
"""
import unittest
import pyopenms


class TestFileTypesStaticMethods(unittest.TestCase):
    """Test that static methods work correctly without instance"""

    def test_typeToName_is_static(self):
        """Test typeToName can be called as static method"""
        # Should be callable without creating an instance
        name = pyopenms.FileTypes.typeToName(pyopenms.FileType.MZML)
        self.assertIsInstance(name, str)
        self.assertEqual(name, "mzML")

    def test_typeToDescription_is_static(self):
        """Test typeToDescription can be called as static method"""
        # Should be callable without creating an instance
        desc = pyopenms.FileTypes.typeToDescription(pyopenms.FileType.MZML)
        self.assertIsInstance(desc, str)
        self.assertEqual(desc, "mzML raw data file")

    def test_nameToType_is_static(self):
        """Test nameToType can be called as static method"""
        # Should be callable without creating an instance
        file_type = pyopenms.FileTypes.nameToType("mzML")
        self.assertEqual(file_type, pyopenms.FileType.MZML)

    def test_typeToMZML_is_static(self):
        """Test typeToMZML can be called as static method"""
        # Should be callable without creating an instance
        mzml_name = pyopenms.FileTypes.typeToMZML(pyopenms.FileType.MZML)
        self.assertIsInstance(mzml_name, str)
        self.assertEqual(mzml_name, "mzML file")


class TestFileTypesEnum(unittest.TestCase):
    """Test FileType enum values"""

    def test_all_file_types_are_accessible(self):
        """Test that all FileType enum values are accessible"""
        file_types = [
            pyopenms.FileType.UNKNOWN,
            pyopenms.FileType.DTA,
            pyopenms.FileType.DTA2D,
            pyopenms.FileType.MZDATA,
            pyopenms.FileType.MZXML,
            pyopenms.FileType.FEATUREXML,
            pyopenms.FileType.IDXML,
            pyopenms.FileType.CONSENSUSXML,
            pyopenms.FileType.MGF,
            pyopenms.FileType.INI,
            pyopenms.FileType.TOPPAS,
            pyopenms.FileType.TRANSFORMATIONXML,
            pyopenms.FileType.MZML,
            pyopenms.FileType.CACHEDMZML,
            pyopenms.FileType.MS2,
            pyopenms.FileType.PEPXML,
            pyopenms.FileType.PROTXML,
            pyopenms.FileType.MZIDENTML,
            pyopenms.FileType.QCML,
            pyopenms.FileType.MZQC,
            pyopenms.FileType.GELML,
            pyopenms.FileType.TRAML,
            pyopenms.FileType.MSP,
            pyopenms.FileType.OMSSAXML,
            pyopenms.FileType.MASCOTXML,
            pyopenms.FileType.PNG,
            pyopenms.FileType.XMASS,
            pyopenms.FileType.TSV,
            pyopenms.FileType.MZTAB,
            pyopenms.FileType.PEPLIST,
            pyopenms.FileType.HARDKLOER,
            pyopenms.FileType.KROENIK,
            pyopenms.FileType.FASTA,
            pyopenms.FileType.EDTA,
            pyopenms.FileType.CSV,
            pyopenms.FileType.TXT,
            pyopenms.FileType.OBO,
            pyopenms.FileType.HTML,
            pyopenms.FileType.ANALYSISXML,
            pyopenms.FileType.XSD,
            pyopenms.FileType.PSQ,
            pyopenms.FileType.MRM,
            pyopenms.FileType.SQMASS,
            pyopenms.FileType.PQP,
            pyopenms.FileType.MS,
            pyopenms.FileType.OSW,
            pyopenms.FileType.PSMS,
            pyopenms.FileType.PIN,
            pyopenms.FileType.PARAMXML,
            pyopenms.FileType.SPLIB,
            pyopenms.FileType.NOVOR,
            pyopenms.FileType.XQUESTXML,
            pyopenms.FileType.SPECXML,
            pyopenms.FileType.JSON,
            pyopenms.FileType.RAW,
            pyopenms.FileType.OMS,
            pyopenms.FileType.EXE,
            pyopenms.FileType.XML,
            pyopenms.FileType.BZ2,
            pyopenms.FileType.GZ,
            pyopenms.FileType.PARQUET,
            pyopenms.FileType.SIZE_OF_TYPE,
        ]
        
        # All enum values should be integers
        for ft in file_types:
            self.assertIsInstance(ft, int)

    def test_file_type_enum_values_are_unique(self):
        """Test that FileType enum values are unique (except duplicates are OK)"""
        # Just check a few key values to ensure they're different
        self.assertNotEqual(pyopenms.FileType.MZML, pyopenms.FileType.MZXML)
        self.assertNotEqual(pyopenms.FileType.FEATUREXML, pyopenms.FileType.IDXML)
        self.assertNotEqual(pyopenms.FileType.UNKNOWN, pyopenms.FileType.MZML)


class TestFileTypesConversions(unittest.TestCase):
    """Test conversion methods between names and types"""

    def test_typeToName_returns_correct_names(self):
        """Test typeToName returns correct file extension names"""
        test_cases = [
            (pyopenms.FileType.MZML, "mzML"),
            (pyopenms.FileType.MZXML, "mzXML"),
            (pyopenms.FileType.FEATUREXML, "featureXML"),
            (pyopenms.FileType.IDXML, "idXML"),
            (pyopenms.FileType.FASTA, "fasta"),
            (pyopenms.FileType.TSV, "tsv"),
            (pyopenms.FileType.CSV, "csv"),
            (pyopenms.FileType.JSON, "json"),
            (pyopenms.FileType.PARQUET, "parquet"),
        ]
        
        for file_type, expected_name in test_cases:
            with self.subTest(file_type=file_type):
                name = pyopenms.FileTypes.typeToName(file_type)
                self.assertIsInstance(name, str)
                self.assertEqual(name, expected_name)

    def test_nameToType_returns_correct_types(self):
        """Test nameToType converts names to correct types"""
        test_cases = [
            ("mzML", pyopenms.FileType.MZML),
            ("mzXML", pyopenms.FileType.MZXML),
            ("featureXML", pyopenms.FileType.FEATUREXML),
            ("idXML", pyopenms.FileType.IDXML),
            ("fasta", pyopenms.FileType.FASTA),
            ("FASTA", pyopenms.FileType.FASTA),  # Case insensitive
            ("json", pyopenms.FileType.JSON),
            ("parquet", pyopenms.FileType.PARQUET),
            ("pqt", pyopenms.FileType.PARQUET),  # Alternative extension
        ]
        
        for name, expected_type in test_cases:
            with self.subTest(name=name):
                file_type = pyopenms.FileTypes.nameToType(name)
                self.assertEqual(file_type, expected_type)

    def test_nameToType_case_insensitive(self):
        """Test that nameToType is case insensitive"""
        # Test various case combinations
        self.assertEqual(
            pyopenms.FileTypes.nameToType("mzml"),
            pyopenms.FileType.MZML
        )
        self.assertEqual(
            pyopenms.FileTypes.nameToType("MZML"),
            pyopenms.FileType.MZML
        )
        self.assertEqual(
            pyopenms.FileTypes.nameToType("MzMl"),
            pyopenms.FileType.MZML
        )

    def test_nameToType_unknown_returns_unknown(self):
        """Test that unknown file type names return UNKNOWN"""
        unknown_type = pyopenms.FileTypes.nameToType("unknown_extension")
        self.assertEqual(unknown_type, pyopenms.FileType.UNKNOWN)

    def test_typeToDescription_returns_descriptions(self):
        """Test typeToDescription returns human-readable descriptions"""
        # Test a few key descriptions
        desc = pyopenms.FileTypes.typeToDescription(pyopenms.FileType.MZML)
        self.assertIsInstance(desc, str)
        self.assertIn("mzML", desc)
        
        desc = pyopenms.FileTypes.typeToDescription(pyopenms.FileType.FEATUREXML)
        self.assertIn("OpenMS", desc)

    def test_typeToMZML_returns_mzml_names(self):
        """Test typeToMZML returns appropriate mzML CV term names"""
        # Only certain types have mzML names
        test_cases = [
            (pyopenms.FileType.MZML, "mzML file"),
            (pyopenms.FileType.MZDATA, "PSI mzData file"),
            (pyopenms.FileType.MZXML, "ISB mzXML file"),
            (pyopenms.FileType.MGF, "Mascot MGF file"),
        ]
        
        for file_type, expected_name in test_cases:
            with self.subTest(file_type=file_type):
                name = pyopenms.FileTypes.typeToMZML(file_type)
                self.assertIsInstance(name, str)
                self.assertEqual(name, expected_name)


class TestFilePropertiesEnum(unittest.TestCase):
    """Test FileProperties enum class"""

    def test_file_properties_enum_accessible(self):
        """Test that FileProperties enum class values are accessible"""
        properties = [
            pyopenms.FileProperties.READABLE,
            pyopenms.FileProperties.WRITEABLE,
            pyopenms.FileProperties.PROVIDES_SPECTRUM,
            pyopenms.FileProperties.PROVIDES_EXPERIMENT,
            pyopenms.FileProperties.PROVIDES_FEATURES,
            pyopenms.FileProperties.PROVIDES_CONSENSUSFEATURES,
            pyopenms.FileProperties.PROVIDES_IDENTIFICATIONS,
            pyopenms.FileProperties.PROVIDES_TRANSITIONS,
            pyopenms.FileProperties.PROVIDES_QUANTIFICATIONS,
            pyopenms.FileProperties.PROVIDES_TRANSFORMATIONS,
            pyopenms.FileProperties.PROVIDES_QC,
            pyopenms.FileProperties.SIZE_OF_FILEPROPERTIES,
        ]
        
        # All enum values should be integers
        for prop in properties:
            self.assertIsInstance(prop, int)

    def test_file_properties_values_unique(self):
        """Test that FileProperties enum values are unique"""
        properties = [
            pyopenms.FileProperties.READABLE,
            pyopenms.FileProperties.WRITEABLE,
            pyopenms.FileProperties.PROVIDES_SPECTRUM,
            pyopenms.FileProperties.PROVIDES_EXPERIMENT,
            pyopenms.FileProperties.PROVIDES_FEATURES,
        ]
        
        # Check that values are unique
        self.assertEqual(len(properties), len(set(properties)))


class TestFilterLayoutEnum(unittest.TestCase):
    """Test FilterLayout enum class"""

    def test_filter_layout_enum_accessible(self):
        """Test that FilterLayout enum class values are accessible"""
        layouts = [
            pyopenms.FilterLayout.COMPACT,
            pyopenms.FilterLayout.ONE_BY_ONE,
            pyopenms.FilterLayout.BOTH,
        ]
        
        # All enum values should be integers
        for layout in layouts:
            self.assertIsInstance(layout, int)

    def test_filter_layout_values_unique(self):
        """Test that FilterLayout enum values are unique"""
        self.assertNotEqual(
            pyopenms.FilterLayout.COMPACT,
            pyopenms.FilterLayout.ONE_BY_ONE
        )
        self.assertNotEqual(
            pyopenms.FilterLayout.COMPACT,
            pyopenms.FilterLayout.BOTH
        )
        self.assertNotEqual(
            pyopenms.FilterLayout.ONE_BY_ONE,
            pyopenms.FilterLayout.BOTH
        )


class TestFileTypesRoundTrip(unittest.TestCase):
    """Test round-trip conversions between types and names"""

    def test_type_to_name_to_type_roundtrip(self):
        """Test that converting type -> name -> type gives back original"""
        test_types = [
            pyopenms.FileType.MZML,
            pyopenms.FileType.MZXML,
            pyopenms.FileType.FEATUREXML,
            pyopenms.FileType.IDXML,
            pyopenms.FileType.FASTA,
            pyopenms.FileType.JSON,
            pyopenms.FileType.PARQUET,
        ]
        
        for original_type in test_types:
            with self.subTest(file_type=original_type):
                name = pyopenms.FileTypes.typeToName(original_type)
                recovered_type = pyopenms.FileTypes.nameToType(name)
                self.assertEqual(recovered_type, original_type)


class TestFileTypesStringInput(unittest.TestCase):
    """Test that nameToType accepts both str and bytes (backward compatibility)"""

    def test_nameToType_with_str(self):
        """Test nameToType accepts str input"""
        file_type = pyopenms.FileTypes.nameToType("mzML")
        self.assertEqual(file_type, pyopenms.FileType.MZML)

    def test_nameToType_with_bytes(self):
        """Test nameToType accepts bytes input (backward compatible)"""
        file_type = pyopenms.FileTypes.nameToType(b"mzML")
        self.assertEqual(file_type, pyopenms.FileType.MZML)

    def test_typeToName_returns_str(self):
        """Test typeToName returns str (not bytes)"""
        name = pyopenms.FileTypes.typeToName(pyopenms.FileType.MZML)
        self.assertIsInstance(name, str)
        self.assertEqual(name, "mzML")

    def test_typeToDescription_returns_str(self):
        """Test typeToDescription returns str (not bytes)"""
        desc = pyopenms.FileTypes.typeToDescription(pyopenms.FileType.MZML)
        self.assertIsInstance(desc, str)
        self.assertIn("mzML", desc)

    def test_typeToMZML_returns_str(self):
        """Test typeToMZML returns str (not bytes)"""
        mzml_name = pyopenms.FileTypes.typeToMZML(pyopenms.FileType.MZML)
        self.assertIsInstance(mzml_name, str)
        self.assertEqual(mzml_name, "mzML file")


if __name__ == "__main__":
    unittest.main()
