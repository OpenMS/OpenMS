import pyopenms

def testSpectrumNativeIDParser():
    """Test the SpectrumNativeIDParser class for parsing spectrum native IDs."""

    # Test isNativeID - recognized native ID prefixes
    assert pyopenms.SpectrumNativeIDParser.isNativeID("scan=123")
    assert pyopenms.SpectrumNativeIDParser.isNativeID("scanId=456")
    assert pyopenms.SpectrumNativeIDParser.isNativeID("scanID=456")  # both cases supported
    assert pyopenms.SpectrumNativeIDParser.isNativeID("controllerType=0 controllerNumber=1 scan=100")
    assert pyopenms.SpectrumNativeIDParser.isNativeID("function=2 process=1 scan=100")
    assert pyopenms.SpectrumNativeIDParser.isNativeID("sample=1 period=1 cycle=42 experiment=1")
    assert pyopenms.SpectrumNativeIDParser.isNativeID("index=789")
    assert pyopenms.SpectrumNativeIDParser.isNativeID("spectrum=101112")
    assert pyopenms.SpectrumNativeIDParser.isNativeID("file=42")

    # Test isNativeID - non-native IDs
    assert not pyopenms.SpectrumNativeIDParser.isNativeID("123")
    assert not pyopenms.SpectrumNativeIDParser.isNativeID("")
    assert not pyopenms.SpectrumNativeIDParser.isNativeID("some_random_string")
    assert not pyopenms.SpectrumNativeIDParser.isNativeID("SCAN=123")  # case-sensitive

    # Test getRegExFromNativeID
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("scan=123") == "scan=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("controllerType=0 scan=100") == "scan=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("function=2 scan=100") == "scan=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("index=456") == "index=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("scanId=789") == "scanId=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("scanID=789") == "scanID=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("spectrum=101") == "spectrum=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("file=42") == "file=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("123") == "(?<GROUP>\\d+)"

    # Test extractScanNumber with CV accession - Thermo format (MS:1000768)
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=42", "MS:1000768") == 42
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=0", "MS:1000768") == 0
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("controllerType=0 controllerNumber=1 scan=123", "MS:1000768") == 123

    # Test extractScanNumber - Waters format (MS:1000769)
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=42", "MS:1000769") == 42

    # Test extractScanNumber - WIFF format (MS:1000770) - cycle * 1000 + experiment
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("sample=1 period=1 cycle=42 experiment=1", "MS:1000770") == 42001
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("sample=1 period=1 cycle=100 experiment=50", "MS:1000770") == 100050

    # Test extractScanNumber - Bruker formats
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=42", "MS:1000771") == 42  # Bruker BAF
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=42", "MS:1000772") == 42  # Bruker U2
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("file=42", "MS:1000773") == 42  # Bruker FID
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("file=42", "MS:1000775") == 42  # Single peak list
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=42", "MS:1000776") == 42  # Thermo/Bruker TDF

    # Test extractScanNumber - index format (MS:1000774) - returns index+1
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("index=42", "MS:1000774") == 43
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("index=0", "MS:1000774") == 1

    # Test extractScanNumber - spectrum format (MS:1000777)
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("spectrum=42", "MS:1000777") == 42

    # Test extractScanNumber - Agilent MassHunter format (MS:1001508)
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scanId=42", "MS:1001508") == 42

    # Test extractScanNumber - mzML unique identifier (MS:1001530)
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("42", "MS:1001530") == 42

    # Test extractScanNumber - invalid accession returns -1
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=42", "MS:9999999") == -1
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=42", "invalid") == -1
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("scan=42", "") == -1

    # Test extractScanNumber - invalid native_id for given accession returns -1
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("spectrum=42", "MS:1000768") == -1
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("", "MS:1000768") == -1
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("malformed", "MS:1000768") == -1

    # Test merged spectra - should return last scan number
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber(
        "controllerType=0 controllerNumber=1 scan=100 merged controllerType=0 controllerNumber=1 scan=200",
        "MS:1000768"
    ) == 200

    # Bruker timsTOF .d (TDF) native IDs use "frame=<F> scan=<S>" with optional trailing
    # tokens (precursor=, windowGroup=, scanStart=/scanEnd=, merged=) emitted by the OpenMS
    # DDA/DIA-PASEF reader (MS:1002818). The parser recognizes the frame= prefix and targets
    # the scan=<int> token as the scan-number proxy. The accession-based cases above stop at
    # scan=/MS:1000776 and never exercise frame=; these pin that contract on the Python side,
    # mirroring the C++ SpectrumNativeIDParser_test "[EXTRA] Bruker TDF frame= / MS:1002818" section.

    # isNativeID recognizes the frame= prefix
    assert pyopenms.SpectrumNativeIDParser.isNativeID("frame=2 scan=529")
    assert pyopenms.SpectrumNativeIDParser.isNativeID("frame=10 scan=42 precursor=3")

    # getRegExFromNativeID targets the scan= token for the frame= format
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("frame=2 scan=529") == "scan=(?<GROUP>\\d+)"
    assert pyopenms.SpectrumNativeIDParser.getRegExFromNativeID("frame=10 scan=42 precursor=3") == "scan=(?<GROUP>\\d+)"

    # extractScanNumber via the Bruker TDF accession (MS:1002818) returns the scan= integer
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("frame=2 scan=529", "MS:1002818") == 529
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("frame=10 scan=42 precursor=3", "MS:1002818") == 42

    # merged PASEF spectra carry multiple scan= tokens; the last one is the merged scan number
    assert pyopenms.SpectrumNativeIDParser.extractScanNumber("frame=1 scan=5 frame=1 scan=9", "MS:1002818") == 9

    print("All SpectrumNativeIDParser tests passed!")


if __name__ == "__main__":
    testSpectrumNativeIDParser()
