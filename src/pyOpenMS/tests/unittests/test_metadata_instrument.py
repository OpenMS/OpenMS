#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Metadata and Instrument tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms

# Import shared test helper functions
from test_helpers import report, _testMetaInfoInterface


@report
def testAcquisitionInfo():
    """
    @tests: AcquisitionInfo
     AcquisitionInfo.__init__
     AcquisitionInfo.__eq__
     AcquisitionInfo.__ge__
     AcquisitionInfo.__gt__
     AcquisitionInfo.__le__
     AcquisitionInfo.__lt__
     AcquisitionInfo.__ne__
     AcquisitionInfo.getMethodOfCombination
     AcquisitionInfo.setMethodOfCombination
    """

    ai = pyopenms.AcquisitionInfo()
    ai.__doc__

    assert ai == ai
    assert not ai != ai
    ai.setMethodOfCombination("ABC")
    assert ai.getMethodOfCombination() == "ABC"


@report
def testInstrumentSettings():
    """
    @tests: InstrumentSettings
     InstrumentSettings.__init__
     InstrumentSettings.clearMetaInfo
     InstrumentSettings.getKeys
     InstrumentSettings.getMetaValue
     InstrumentSettings.getPolarity
     InstrumentSettings.isMetaEmpty
     InstrumentSettings.metaValueExists
     InstrumentSettings.removeMetaValue
     InstrumentSettings.setMetaValue
     InstrumentSettings.setPolarity
     InstrumentSettings.__eq__
     InstrumentSettings.__ge__
     InstrumentSettings.__gt__
     InstrumentSettings.__le__
     InstrumentSettings.__lt__
     InstrumentSettings.__ne__
     """
    ins = pyopenms.InstrumentSettings()
    _testMetaInfoInterface(ins)
    ins.setPolarity(pyopenms.IonSource.Polarity.NEGATIVE)
    assert ins.getPolarity() == pyopenms.IonSource.Polarity.NEGATIVE

    assert ins == ins
    assert not ins != ins

    # Test getAllNamesOf method
    scan_mode_names = pyopenms.InstrumentSettings.getAllNamesOfScanMode()
    assert len(scan_mode_names) == pyopenms.InstrumentSettings.ScanMode.SIZE_OF_SCANMODE
    assert scan_mode_names[pyopenms.InstrumentSettings.ScanMode.MS1SPECTRUM].decode() == "MS1Spectrum"


@report
def testContactPerson():
    """
    @tests: ContactPerson
     ContactPerson.__init__
     ContactPerson.getFirstName
     ContactPerson.setFirstName
     ContactPerson.getLastName
     ContactPerson.setLastName
     ContactPerson.setName
     ContactPerson.getInstitution
     ContactPerson.setInstitution
     ContactPerson.getEmail
     ContactPerson.setEmail
     ContactPerson.getURL
     ContactPerson.setURL
     ContactPerson.getAddress
     ContactPerson.setAddress
     ContactPerson.getContactInfo
     ContactPerson.setContactInfo
     """
    ins = pyopenms.ContactPerson()

    ins.getFirstName()
    ins.setFirstName("test")
    ins.getLastName()
    ins.setLastName("test")
    ins.setName("Testy Test")
    ins.getInstitution()
    ins.setInstitution("test")
    ins.getEmail()
    ins.setEmail("test")
    ins.getURL()
    ins.setURL("test")
    ins.getAddress()
    ins.setAddress("test")
    ins.getContactInfo()
    ins.setContactInfo("test")


@report
def testDocumentIdentifier():
    """
    @tests: DocumentIdentifier
     DocumentIdentifier.__init__
     DocumentIdentifier.setIdentifier
     DocumentIdentifier.getIdentifier
     DocumentIdentifier.setLoadedFilePath
     DocumentIdentifier.getLoadedFilePath
     DocumentIdentifier.setLoadedFileType
     DocumentIdentifier.getLoadedFileType
     """
    ins = pyopenms.DocumentIdentifier()

    ins.setIdentifier("test")
    ins.getIdentifier()
    # ins.setLoadedFilePath("Test")
    ins.getLoadedFilePath()
    # ins.setLoadedFileType("test")
    ins.getLoadedFileType()


@report
def testGradient():
    """
    @tests: Gradient
     Gradient.__init__
     Gradient.addEluent
     Gradient.addEluent
     Gradient.clearEluents
     Gradient.getEluents
     Gradient.addTimepoint
     Gradient.clearTimepoints
     Gradient.getTimepoints
     Gradient.getPercentage
     Gradient.setPercentage
     Gradient.clearPercentages
     Gradient.isValid
     """
    ins = pyopenms.Gradient()

    ins.addEluent("test")
    ins.clearEluents()
    assert len(ins.getEluents() ) == 0
    ins.addEluent("test")
    assert len(ins.getEluents() ) == 1

    ins.clearTimepoints()
    ins.addTimepoint(5)
    ins.getTimepoints()

    ins.setPercentage("test", 5, 20)
    ins.getPercentage("test", 5)
    ins.clearPercentages()
    ins.isValid()


@report
def testHPLC():
    """
    @tests: HPLC
     HPLC.__init__
     HPLC.getInstrument
     HPLC.setInstrument
     HPLC.getColumn
     HPLC.setColumn
     HPLC.getTemperature
     HPLC.setTemperature
     HPLC.getPressure
     HPLC.setPressure
     HPLC.getFlux
     HPLC.setFlux
     HPLC.setComment
     HPLC.getComment
     HPLC.setGradient
     HPLC.getGradient
     """
    ins = pyopenms.HPLC()

    ins.setInstrument("test")
    ins.getInstrument()
    ins.setColumn("test")
    ins.getColumn()
    ins.setTemperature(6)
    ins.getTemperature()
    ins.setPressure(6)
    ins.getPressure()
    ins.setFlux(8)
    ins.getFlux()
    ins.setComment("test")
    ins.getComment()

    g = pyopenms.Gradient()
    ins.setGradient(g)
    ins.getGradient()


@report
def testInstrument():
    """
    @tests: Instrument
     Instrument.__init__
     Instrument.setName
     Instrument.getName
     Instrument.setVendor
     Instrument.getVendor
     Instrument.setModel
     Instrument.getModel
     Instrument.setCustomizations
     Instrument.getCustomizations
     Instrument.setIonSources
     Instrument.getIonSources
     Instrument.setMassAnalyzers
     Instrument.getMassAnalyzers
     Instrument.setIonDetectors
     Instrument.getIonDetectors
     Instrument.setSoftware
     Instrument.getSoftware
     """
    ins = pyopenms.Instrument()

    ins.setName("test")
    ins.getName()
    ins.setVendor("test")
    ins.getVendor()
    ins.setModel("test")
    ins.getModel()
    ins.setCustomizations("test")
    ins.getCustomizations()

    ion_sources = [ pyopenms.IonSource() for i in range(5)]
    ins.setIonSources(ion_sources)
    ins.getIonSources()
    mass_analyzers = [ pyopenms.MassAnalyzer() for i in range(5)]
    ins.setMassAnalyzers(mass_analyzers)
    ins.getMassAnalyzers()
    ion_detectors = [ pyopenms.IonDetector() for i in range(5)]
    ins.setIonDetectors(ion_detectors)
    ins.getIonDetectors()

    s = pyopenms.Software()
    ins.setSoftware(s)
    ins.getSoftware()

    # Test getAllNamesOf method
    ion_optics_names = pyopenms.Instrument.getAllNamesOfIonOpticsType()
    assert len(ion_optics_names) == pyopenms.Instrument.IonOpticsType.SIZE_OF_IONOPTICSTYPE
    assert ion_optics_names[pyopenms.Instrument.IonOpticsType.REFLECTRON].decode() == "reflectron"


@report
def testIonDetector():
    """
    @tests: IonDetector
     IonDetector.__init__
     IonDetector.setAcquisitionMode
     IonDetector.getAcquisitionMode
     IonDetector.setResolution
     IonDetector.getResolution
     IonDetector.setADCSamplingFrequency
     IonDetector.getADCSamplingFrequency
     IonDetector.setOrder
     IonDetector.getOrder
     IonDetector.getAllNamesOfType
     IonDetector.getAllNamesOfAcquisitionMode
     """
    ins = pyopenms.IonDetector()

    m = pyopenms.IonDetector.AcquisitionMode.ACQMODENULL
    ins.setAcquisitionMode(m)
    ins.getAcquisitionMode()

    ins.setResolution(8.0)
    ins.getResolution()

    ins.setADCSamplingFrequency(8.0)
    ins.getADCSamplingFrequency()

    ins.setOrder(8)
    ins.getOrder()

    # Test getAllNamesOf methods
    type_names = pyopenms.IonDetector.getAllNamesOfType()
    assert len(type_names) == pyopenms.IonDetector.Type_IonDetector.SIZE_OF_TYPE
    assert type_names[pyopenms.IonDetector.Type_IonDetector.ELECTRONMULTIPLIER].decode() == "Electron multiplier"

    acq_mode_names = pyopenms.IonDetector.getAllNamesOfAcquisitionMode()
    assert len(acq_mode_names) == pyopenms.IonDetector.AcquisitionMode.SIZE_OF_ACQUISITIONMODE
    assert acq_mode_names[pyopenms.IonDetector.AcquisitionMode.PULSECOUNTING].decode() == "Pulse counting"


@report
def testIonSource():
    """
    @tests: IonSource
     IonSource.__init__
     IonSource.setPolarity
     IonSource.getPolarity
     IonSource.setInletType
     IonSource.getInletType
     IonSource.setIonizationMethod
     IonSource.getIonizationMethod
     IonSource.setOrder
     IonSource.getOrder
     IonSource.getAllNamesOfInletType
     IonSource.getAllNamesOfIonizationMethod
     IonSource.getAllNamesOfPolarity
     """
    ins = pyopenms.IonSource()

    p = pyopenms.IonSource.Polarity.POSITIVE
    ins.setPolarity(p)
    ins.getPolarity()

    i = pyopenms.IonSource.InletType.INLETNULL
    ins.setInletType(i)
    ins.getInletType()

    i = pyopenms.IonSource.IonizationMethod.ESI
    ins.setIonizationMethod(i)
    ins.getIonizationMethod()

    ins.setOrder(5)
    ins.getOrder()

    # Test getAllNamesOf methods
    inlet_names = pyopenms.IonSource.getAllNamesOfInletType()
    assert len(inlet_names) == pyopenms.IonSource.InletType.SIZE_OF_INLETTYPE
    assert inlet_names[pyopenms.IonSource.InletType.INLETNULL].decode() == "Unknown"
    assert inlet_names[pyopenms.IonSource.InletType.DIRECT].decode() == "Direct"
    assert inlet_names[pyopenms.IonSource.InletType.NANOSPRAY].decode() == "Nanospray inlet"

    ionization_names = pyopenms.IonSource.getAllNamesOfIonizationMethod()
    assert len(ionization_names) == pyopenms.IonSource.IonizationMethod.SIZE_OF_IONIZATIONMETHOD
    assert ionization_names[pyopenms.IonSource.IonizationMethod.IONMETHODNULL].decode() == "Unknown"
    assert ionization_names[pyopenms.IonSource.IonizationMethod.ESI].decode() == "Electrospray ionisation"
    assert ionization_names[pyopenms.IonSource.IonizationMethod.MALDI].decode() == "Matrix-assisted laser desorption ionization"

    polarity_names = pyopenms.IonSource.getAllNamesOfPolarity()
    assert len(polarity_names) == pyopenms.IonSource.Polarity.SIZE_OF_POLARITY
    assert polarity_names[pyopenms.IonSource.Polarity.POLNULL].decode() == "unknown"
    assert polarity_names[pyopenms.IonSource.Polarity.POSITIVE].decode() == "positive"
    assert polarity_names[pyopenms.IonSource.Polarity.NEGATIVE].decode() == "negative"


@report
def testMassAnalyzer():
    """
    @tests: MassAnalyzer
     MassAnalyzer.__init__
     MassAnalyzer.setType
     MassAnalyzer.getType
     MassAnalyzer.setResolutionMethod
     MassAnalyzer.getResolutionMethod
     MassAnalyzer.setResolutionType
     MassAnalyzer.getResolutionType
     MassAnalyzer.setScanDirection
     MassAnalyzer.getScanDirection
     MassAnalyzer.setScanLaw
     MassAnalyzer.getScanLaw
     MassAnalyzer.setReflectronState
     MassAnalyzer.getReflectronState
     MassAnalyzer.setResolution
     MassAnalyzer.getResolution
     MassAnalyzer.setAccuracy
     MassAnalyzer.getAccuracy
     MassAnalyzer.setScanRate
     MassAnalyzer.getScanRate
     MassAnalyzer.setScanTime
     MassAnalyzer.getScanTime
     MassAnalyzer.setTOFTotalPathLength
     MassAnalyzer.getTOFTotalPathLength
     MassAnalyzer.setIsolationWidth
     MassAnalyzer.getIsolationWidth
     MassAnalyzer.setFinalMSExponent
     MassAnalyzer.getFinalMSExponent
     MassAnalyzer.setMagneticFieldStrength
     MassAnalyzer.getMagneticFieldStrength
     MassAnalyzer.setOrder
     MassAnalyzer.getOrder
     MassAnalyzer.getAllNamesOfAnalyzerType
     MassAnalyzer.getAllNamesOfResolutionMethod
     MassAnalyzer.getAllNamesOfResolutionType
     MassAnalyzer.getAllNamesOfScanDirection
     MassAnalyzer.getAllNamesOfScanLaw
     MassAnalyzer.getAllNamesOfReflectronState
     """
    ins = pyopenms.MassAnalyzer()

    ma = pyopenms.MassAnalyzer.AnalyzerType.QUADRUPOLE
    ins.setType(ma)
    ins.getType()

    res = pyopenms.MassAnalyzer.ResolutionMethod.FWHM
    ins.setResolutionMethod(res)
    ins.getResolutionMethod()

    res = pyopenms.MassAnalyzer.ResolutionType.CONSTANT
    ins.setResolutionType(res)
    ins.getResolutionType()

    res = pyopenms.MassAnalyzer.ScanDirection.UP
    ins.setScanDirection(res)
    ins.getScanDirection()

    res = pyopenms.MassAnalyzer.ScanLaw.LINEAR
    ins.setScanLaw(res)
    ins.getScanLaw()

    res = pyopenms.MassAnalyzer.ReflectronState.ON
    ins.setReflectronState(res)
    ins.getReflectronState()

    ins.setResolution(5.0)
    ins.getResolution()
    ins.setAccuracy(5.0)
    ins.getAccuracy()
    ins.setScanRate(5.0)
    ins.getScanRate()
    ins.setScanTime(5.0)
    ins.getScanTime()
    ins.setTOFTotalPathLength(5.0)
    ins.getTOFTotalPathLength()
    ins.setIsolationWidth(5.0)
    ins.getIsolationWidth()
    ins.setFinalMSExponent(5)
    ins.getFinalMSExponent()
    ins.setMagneticFieldStrength(5.0)
    ins.getMagneticFieldStrength()
    ins.setOrder(5)
    ins.getOrder()

    # Test getAllNamesOf methods
    analyzer_names = pyopenms.MassAnalyzer.getAllNamesOfAnalyzerType()
    assert len(analyzer_names) == pyopenms.MassAnalyzer.AnalyzerType.SIZE_OF_ANALYZERTYPE
    assert analyzer_names[pyopenms.MassAnalyzer.AnalyzerType.QUADRUPOLE].decode() == "Quadrupole"
    assert analyzer_names[pyopenms.MassAnalyzer.AnalyzerType.ORBITRAP].decode() == "Orbitrap"

    res_method_names = pyopenms.MassAnalyzer.getAllNamesOfResolutionMethod()
    assert len(res_method_names) == pyopenms.MassAnalyzer.ResolutionMethod.SIZE_OF_RESOLUTIONMETHOD
    assert res_method_names[pyopenms.MassAnalyzer.ResolutionMethod.FWHM].decode() == "Full width at half max"

    res_type_names = pyopenms.MassAnalyzer.getAllNamesOfResolutionType()
    assert len(res_type_names) == pyopenms.MassAnalyzer.ResolutionType.SIZE_OF_RESOLUTIONTYPE
    assert res_type_names[pyopenms.MassAnalyzer.ResolutionType.CONSTANT].decode() == "Constant"

    scan_dir_names = pyopenms.MassAnalyzer.getAllNamesOfScanDirection()
    assert len(scan_dir_names) == pyopenms.MassAnalyzer.ScanDirection.SIZE_OF_SCANDIRECTION
    assert scan_dir_names[pyopenms.MassAnalyzer.ScanDirection.UP].decode() == "Up"

    scan_law_names = pyopenms.MassAnalyzer.getAllNamesOfScanLaw()
    assert len(scan_law_names) == pyopenms.MassAnalyzer.ScanLaw.SIZE_OF_SCANLAW
    assert scan_law_names[pyopenms.MassAnalyzer.ScanLaw.LINEAR].decode() == "Linar"  # Note: typo in source

    reflectron_names = pyopenms.MassAnalyzer.getAllNamesOfReflectronState()
    assert len(reflectron_names) == pyopenms.MassAnalyzer.ReflectronState.SIZE_OF_REFLECTRONSTATE
    assert reflectron_names[pyopenms.MassAnalyzer.ReflectronState.ON].decode() == "On"


@report
def testSample():
    """
    @tests: Sample
     Sample.__init__
     Sample.setName
     Sample.getName
     Sample.setOrganism
     Sample.getOrganism
     Sample.setNumber
     Sample.getNumber
     Sample.setComment
     Sample.getComment
     Sample.setState
     Sample.getState
     Sample.setMass
     Sample.getMass
     Sample.setVolume
     Sample.getVolume
     Sample.setConcentration
     Sample.getConcentration
     Sample.getSubsamples
     Sample.setSubsamples
     """
    ins = pyopenms.Sample()

    ins.setName("test")
    ins.getName()
    ins.setOrganism("test")
    ins.getOrganism()
    ins.setNumber("test")
    ins.getNumber()
    ins.setComment("test")
    ins.getComment()

    state = pyopenms.Sample.SampleState.LIQUID
    ins.setState(state)
    ins.getState()
    ins.setMass(42.0)
    ins.getMass()
    ins.setVolume(42.0)
    ins.getVolume()
    ins.setConcentration(42.0)
    ins.getConcentration()

    a = ins.getSubsamples()
    ins.setSubsamples(a)


@report
def testSoftware():
    """
    @tests: Software
     Software.__init__
     Software.getName
     Software.getVersion
     Software.setName
     Software.setVersion
    """
    sw = pyopenms.Software()
    sw.setName("name")
    sw.setVersion("1.0.0")
    assert sw.getName() == "name"
    assert sw.getVersion() == "1.0.0"


@report
def testSourceFile():
    """
    @tests: SourceFile
     SourceFile.__init__
     SourceFile.getChecksum
     SourceFile.getChecksumType
     SourceFile.getFileSize
     SourceFile.getFileType
     SourceFile.getNameOfFile
     SourceFile.getNativeIDType
     SourceFile.getPathToFile
     SourceFile.setChecksum
     SourceFile.setFileSize
     SourceFile.setFileType
     SourceFile.setNameOfFile
     SourceFile.setNativeIDType
     SourceFile.setPathToFile

    """
    sf = pyopenms.SourceFile()
    sf.setNameOfFile("file.txt")
    assert sf.getNameOfFile() == "file.txt"
    sf.setPathToFile("file.txt")
    assert sf.getPathToFile() == "file.txt"
    sf.setFileType(".txt")
    assert sf.getFileType() == ".txt"
    sf.setChecksum("abcde000", pyopenms.ChecksumType.UNKNOWN_CHECKSUM)
    assert sf.getChecksum() == "abcde000"

    assert sf.getChecksumType() in (pyopenms.ChecksumType.UNKNOWN_CHECKSUM,
                                    pyopenms.ChecksumType.SHA1,
                                    pyopenms.ChecksumType.MD5)

    # Test getAllNamesOf method
    checksum_names = pyopenms.SourceFile.getAllNamesOfChecksumType()
    assert len(checksum_names) == pyopenms.ChecksumType.SIZE_OF_CHECKSUMTYPE
    assert checksum_names[pyopenms.ChecksumType.SHA1].decode() == "SHA-1"


@report
def testDataProcessing():
    """
    @tests: DataProcessing
     DataProcessing.__init__
     DataProcessing.getKeys
     DataProcessing.getMetaValue
     DataProcessing.getProcessingActions
     DataProcessing.getSoftware
     DataProcessing.isMetaEmpty
     DataProcessing.metaValueExists
     DataProcessing.removeMetaValue
     DataProcessing.setCompletionTime
     DataProcessing.setMetaValue
     DataProcessing.setProcessingActions
     DataProcessing.setSoftware
     DataProcessing.__eq__
     DataProcessing.__ge__
     DataProcessing.__gt__
     DataProcessing.__le__
     DataProcessing.__lt__
     DataProcessing.__ne__
     DataProcessing.clearMetaInfo
     DataProcessing.getCompletionTime
    """
    dp = pyopenms.DataProcessing()

    _testMetaInfoInterface(dp)

    assert dp == dp
    assert not dp != dp

    # assert isinstance(dp.getCompletionTime().getDate(), bytes)
    # assert isinstance(dp.getCompletionTime().getTime(), bytes)
    dp.clearMetaInfo()


if __name__ == "__main__":
    # Run all tests in this module
    import sys
    module = sys.modules[__name__]
    
    test_functions = [getattr(module, name) for name in dir(module) 
                     if name.startswith('test') and callable(getattr(module, name))]
    
    print(f"Running {len(test_functions)} metadata/instrument tests...")
    for test_func in test_functions:
        try:
            test_func()
            print(f"✅ {test_func.__name__}")
        except Exception as e:
            print(f"❌ {test_func.__name__}: {e}")
