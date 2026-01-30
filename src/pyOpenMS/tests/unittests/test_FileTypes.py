import pyopenms


def test_enum_values_accessible():
    """Test that all FileType enum values are accessible."""
    expected = [
        "UNKNOWN", "DTA", "DTA2D", "MZDATA", "MZXML", "FEATUREXML", "IDXML",
        "CONSENSUSXML", "MGF", "INI", "TOPPAS", "TRANSFORMATIONXML", "MZML",
        "CACHEDMZML", "MS2", "PEPXML", "PROTXML", "MZIDENTML", "QCML", "MZQC",
        "GELML", "TRAML", "MSP", "OMSSAXML", "MASCOTXML", "PNG", "XMASS",
        "TSV", "MZTAB", "PEPLIST", "HARDKLOER", "KROENIK", "FASTA", "EDTA",
        "CSV", "TXT", "OBO", "HTML", "ANALYSISXML", "XSD", "PSQ", "MRM",
        "SQMASS", "PQP", "MS", "OSW", "PSMS", "PIN", "PARAMXML", "SPLIB",
        "NOVOR", "XQUESTXML", "SPECXML", "JSON", "RAW", "OMS", "EXE", "XML",
        "BZ2", "GZ", "PARQUET", "SIZE_OF_TYPE",
    ]
    for name in expected:
        val = getattr(pyopenms.FileType, name)
        assert val is not None


def test_typeToName():
    """Test typeToName returns correct strings for known types."""
    ft = pyopenms.FileTypes()
    assert ft.typeToName(pyopenms.FileType.MZML) == "mzML"
    assert ft.typeToName(pyopenms.FileType.FASTA) == "fasta"
    assert ft.typeToName(pyopenms.FileType.IDXML) == "idXML"
    assert ft.typeToName(pyopenms.FileType.FEATUREXML) == "featureXML"


def test_typeToDescription():
    """Test typeToDescription returns non-empty strings."""
    ft = pyopenms.FileTypes()
    desc = ft.typeToDescription(pyopenms.FileType.MZML)
    assert isinstance(desc, str)
    assert len(desc) > 0

    desc2 = ft.typeToDescription(pyopenms.FileType.FASTA)
    assert isinstance(desc2, str)
    assert len(desc2) > 0


def test_nameToType():
    """Test nameToType converts names back to types."""
    ft = pyopenms.FileTypes()
    assert ft.nameToType("mzML") == pyopenms.FileType.MZML
    assert ft.nameToType("fasta") == pyopenms.FileType.FASTA
    assert ft.nameToType("idXML") == pyopenms.FileType.IDXML


def test_nameToType_roundtrip():
    """Test that typeToName and nameToType round-trip correctly."""
    ft = pyopenms.FileTypes()
    for t in [pyopenms.FileType.MZML, pyopenms.FileType.FASTA,
              pyopenms.FileType.TRAML, pyopenms.FileType.FEATUREXML]:
        name = ft.typeToName(t)
        assert ft.nameToType(name) == t


def test_typeToMZML():
    """Test typeToMZML for known types."""
    ft = pyopenms.FileTypes()
    result = ft.typeToMZML(pyopenms.FileType.MZML)
    assert isinstance(result, str)


def test_static_methods_callable_on_class():
    """Test that static methods are callable without an instance."""
    # These should work as static methods on the class
    name = pyopenms.FileTypes.typeToName(pyopenms.FileType.MZML)
    assert name == "mzML"

    desc = pyopenms.FileTypes.typeToDescription(pyopenms.FileType.MZML)
    assert len(desc) > 0

    t = pyopenms.FileTypes.nameToType("mzML")
    assert t == pyopenms.FileType.MZML

    mzml_name = pyopenms.FileTypes.typeToMZML(pyopenms.FileType.MZML)
    assert isinstance(mzml_name, str)
