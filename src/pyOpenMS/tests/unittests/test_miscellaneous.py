#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Miscellaneous tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## 
## This file contains various remaining tests that don't fit into other
## categories: data types, enums, RNA sequences, version info, etc.
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms
import copy
import os

# Import shared test helper functions
from test_helpers import report, _testStrOutput


@report
def testConsensusMap():
    """
    @tests: ConsensusMap
     ConsensusMap.__eq__
     ConsensusMap.__ge__
     ConsensusMap.__gt__
     ConsensusMap.__init__
     ConsensusMap.__iter__
     ConsensusMap.__le__
     ConsensusMap.__len__
     ConsensusMap.__lt__
     ConsensusMap.__ne__
     ConsensusMap.append
     ConsensusMap.clear
     ConsensusMap.clearUniqueId
     ConsensusMap.ensureUniqueId
     ConsensusMap.extend
     ConsensusMap.getDataProcessing
     ConsensusMap.getColumnHeaders
     ConsensusMap.getProteinIdentifications
     ConsensusMap.getUnassignedPeptideIdentifications
     ConsensusMap.getUniqueId
     ConsensusMap.hasInvalidUniqueId
     ConsensusMap.hasValidUniqueId
     ConsensusMap.setDataProcessing
     ConsensusMap.setColumnHeaders
     ConsensusMap.setProteinIdentifications
     ConsensusMap.setUnassignedPeptideIdentifications
     ConsensusMap.setUniqueId
     ConsensusMap.setUniqueIds
     ConsensusMap.size
     ConsensusMap.sortByIntensity
     ConsensusMap.sortByMZ
     ConsensusMap.sortByMaps
     ConsensusMap.sortByPosition
     ConsensusMap.sortByQuality
     ConsensusMap.sortByRT
     ConsensusMap.sortBySize
     ConsensusMap.updateRanges
     """
    m = pyopenms.ConsensusMap()
    m_ = copy.copy(m)
    assert m_ == m
    m_ = copy.deepcopy(m)
    assert m_ == m
    m_ = pyopenms.ConsensusMap(m)
    assert m_ == m

    m.clear()
    m.clearUniqueId()
    m.ensureUniqueId()
    m.getDataProcessing()
    m.getColumnHeaders()
    m.getProteinIdentifications()
    m.getUnassignedPeptideIdentifications()
    m.getUniqueId()
    m.hasInvalidUniqueId()
    m.hasValidUniqueId()
    m.setDataProcessing
    m.setColumnHeaders
    m.setProteinIdentifications
    m.setUnassignedPeptideIdentifications
    m.setUniqueId
    m.setUniqueIds
    m.size()
    m.sortByIntensity()
    m.sortByMZ()
    m.sortByMaps()
    m.sortByPosition()
    m.sortByQuality()
    m.sortByRT()
    m.sortBySize()
    # We need to add a feature to the map before calling updateRanges otherwise the getMin and getMax throw an error.
    f = pyopenms.ConsensusFeature()
    f.setCharge(1)
    f.setQuality(2.0)
    f.setWidth(4.0)
    m.push_back(f)
    m.push_back(f)

    m.updateRanges()

    assert isinstance(m.getMinRT(), float)
    assert isinstance(m.getMinRT(), float)
    assert isinstance(m.getMaxMZ(), float)
    assert isinstance(m.getMaxMZ(), float)
    assert isinstance(m.getMinIntensity(), float)
    assert isinstance(m.getMaxIntensity(), float)

    m.getIdentifier()
    m.getLoadedFileType()
    m.getLoadedFilePath()
    
 
    intydf = m.get_intensity_df()
    metadf = m.get_metadata_df()
    assert intydf.shape[0] == 2
    assert metadf.shape[0] == 2

    assert m == m
    assert not m != m

    # Test __repr__ and __str__ methods
    cm_repr = pyopenms.ConsensusMap()
    cf1 = pyopenms.ConsensusFeature()
    cf1.setRT(100.0)
    cf1.setMZ(500.0)
    cf2 = pyopenms.ConsensusFeature()
    cf2.setRT(200.0)
    cf2.setMZ(600.0)
    cm_repr.push_back(cf1)
    cm_repr.push_back(cf2)
    cm_repr.setExperimentType("label-free")

    repr_str = repr(cm_repr)
    assert "ConsensusMap(" in repr_str
    assert "num_consensus_features=" in repr_str

    # Test __len__, append, and extend methods
    cm_len = pyopenms.ConsensusMap()
    assert len(cm_len) == 0
    assert len(cm_len) == cm_len.size()

    cf_test1 = pyopenms.ConsensusFeature()
    cf_test1.setRT(100.0)
    cf_test1.setMZ(500.0)

    cf_test2 = pyopenms.ConsensusFeature()
    cf_test2.setRT(200.0)
    cf_test2.setMZ(600.0)

    cf_test3 = pyopenms.ConsensusFeature()
    cf_test3.setRT(300.0)
    cf_test3.setMZ(700.0)

    # Test append (single item)
    cm_len.append(cf_test1)
    assert len(cm_len) == 1
    assert len(cm_len) == cm_len.size()

    # Test extend with list
    cm_len.extend([cf_test2, cf_test3])
    assert len(cm_len) == 3
    assert len(cm_len) == cm_len.size()

    # Verify the features were added correctly
    assert cm_len[0].getRT() == 100.0
    assert cm_len[1].getRT() == 200.0
    assert cm_len[2].getRT() == 300.0

    # Test extend with another ConsensusMap
    cm_source = pyopenms.ConsensusMap()
    cf_test4 = pyopenms.ConsensusFeature()
    cf_test4.setRT(400.0)
    cf_test4.setMZ(800.0)
    cm_source.push_back(cf_test4)

    cm_len.extend(cm_source)
    assert len(cm_len) == 4
    assert cm_len[3].getRT() == 400.0


@report
def testDataType():
    """
    @tests: DataType
     DataType.DOUBLE_LIST
     DataType.DOUBLE_VALUE
     DataType.EMPTY_VALUE
     DataType.INT_LIST
     DataType.INT_VALUE
     DataType.STRING_LIST
     DataType.STRING_VALUE
    """
    assert isinstance(pyopenms.DataType.DOUBLE_LIST, int)
    assert isinstance(pyopenms.DataType.DOUBLE_VALUE, int)
    assert isinstance(pyopenms.DataType.EMPTY_VALUE, int)
    assert isinstance(pyopenms.DataType.INT_LIST, int)
    assert isinstance(pyopenms.DataType.INT_VALUE, int)
    assert isinstance(pyopenms.DataType.STRING_LIST, int)
    assert isinstance(pyopenms.DataType.STRING_VALUE, int)


@report
def testDataValue():
    """
    @tests: DataValue
     DataValue.__init__
     DataValue.isEmpty
     DataValue.toDoubleList
     DataValue.toDouble
     DataValue.toInt
     DataValue.toIntList
     DataValue.toString
     DataValue.toStringList
     DataValue.valueType

    """
    a = pyopenms.DataValue()
    assert a.isEmpty()

    a = pyopenms.DataValue(1)
    assert not a.isEmpty()
    assert a.toInt() == 1
    assert a.valueType() == pyopenms.DataType.INT_VALUE

    a = pyopenms.DataValue(1.0)
    assert not a.isEmpty()
    assert a.toDouble() == 1.0
    assert a.valueType() == pyopenms.DataType.DOUBLE_VALUE

    a = pyopenms.DataValue("1")
    assert not a.isEmpty()
    assert a.toString() == "1"
    assert a.valueType() == pyopenms.DataType.STRING_VALUE

    a = pyopenms.DataValue([1])
    assert not a.isEmpty()
    assert a.toIntList() == [1]
    assert a.valueType() == pyopenms.DataType.INT_LIST

    a = pyopenms.DataValue([1.0])
    assert not a.isEmpty()
    assert a.toDoubleList() == [1.0]
    assert a.valueType() == pyopenms.DataType.DOUBLE_LIST

    a = pyopenms.DataValue([b"1.0"])
    assert not a.isEmpty()
    assert a.toStringList() == [b"1.0"]
    assert a.valueType() == pyopenms.DataType.STRING_LIST

    assert pyopenms.MSSpectrum().getMetaValue("nonexisingkey") is None


@report
def testLogType():

    """
    @tests: LogType
     LogType.CMD
     LogType.GUI
     LogType.NONE
     """
    assert isinstance(pyopenms.LogType.CMD, int)
    assert isinstance(pyopenms.LogType.GUI, int)
    assert isinstance(pyopenms.LogType.NONE, int)


@report
def testPolarity():
    """
    @tests: Polarity
     Polarity.NEGATIVE
     Polarity.POLNULL
     Polarity.POSITIVE
     Polarity.SIZE_OF_POLARITY
    """
    assert isinstance(pyopenms.IonSource.Polarity.NEGATIVE, int)
    assert isinstance(pyopenms.IonSource.Polarity.POLNULL, int)
    assert isinstance(pyopenms.IonSource.Polarity.POSITIVE, int)


@report
def testPrecursor():
    """
    @tests: Precursor
     Precursor.__init__
     Precursor.getActivationMethods
     Precursor.setActivationMethods
     Precursor.getActivationEnergy
     Precursor.setActivationEnergy
     Precursor.getIsolationWindowLowerOffset
     Precursor.setIsolationWindowLowerOffset
     Precursor.getDriftTime
     Precursor.setDriftTime
     Precursor.getIsolationWindowUpperOffset
     Precursor.setIsolationWindowUpperOffset 
     Precursor.getDriftTimeWindowLowerOffset
     Precursor.setDriftTimeWindowLowerOffset
     Precursor.getDriftTimeWindowUpperOffset
     Precursor.setDriftTimeWindowUpperOffset
     Precursor.getCharge
     Precursor.setCharge
     Precursor.getPossibleChargeStates
     Precursor.setPossibleChargeStates
     Precursor.getUnchargedMass
     Precursor.__eq__
     Precursor.__ne__
    """

    # Test constructors
    prec = pyopenms.Precursor()
    prec2 = pyopenms.Precursor(prec)
    assert prec == prec2
    assert not prec != prec2

    # Test activation methods
    methods = set([pyopenms.Precursor.ActivationMethod.CID, 
                  pyopenms.Precursor.ActivationMethod.HCD])
    methods_short = ["CID", "HCD"]
    methods_long = ["Collision-induced dissociation", "beam-type collision-induced dissociation"]
    prec.setActivationMethods(methods)
    assert prec.getActivationMethods() == methods

    # Test activation energy  
    prec.setActivationEnergy(25.0)
    assert abs(prec.getActivationEnergy() - 25.0) < 1e-5

    # Test activation methods as strings
    short_strings = prec.getActivationMethodsAsShortString()
    assert sorted([s.decode() for s in short_strings]) == sorted(methods_short)
    long_strings = prec.getActivationMethodsAsString()
    assert sorted([s.decode() for s in long_strings]) == sorted(methods_long)

    # Test static methods for all activation methods
    all_names = pyopenms.Precursor.getAllNamesOfActivationMethods()
    assert len(all_names) == pyopenms.Precursor.ActivationMethod.SIZE_OF_ACTIVATIONMETHOD
    assert all_names[pyopenms.Precursor.ActivationMethod.CID].decode() == "Collision-induced dissociation"
    
    all_short_names = pyopenms.Precursor.getAllShortNamesOfActivationMethods()
    assert len(all_short_names) == pyopenms.Precursor.ActivationMethod.SIZE_OF_ACTIVATIONMETHOD
    assert all_short_names[pyopenms.Precursor.ActivationMethod.CID].decode() == "CID"

    # Test isolation window
    prec.setIsolationWindowLowerOffset(0.5)
    assert abs(prec.getIsolationWindowLowerOffset() - 0.5) < 1e-5
    
    prec.setIsolationWindowUpperOffset(1.5) 
    assert abs(prec.getIsolationWindowUpperOffset() - 1.5) < 1e-5

    # Test drift time
    prec.setDriftTime(5.0)
    assert abs(prec.getDriftTime() - 5.0) < 1e-5

    # Test drift time window
    prec.setDriftTimeWindowLowerOffset(0.2)
    assert abs(prec.getDriftTimeWindowLowerOffset() - 0.2) < 1e-5
    
    prec.setDriftTimeWindowUpperOffset(0.8)
    assert abs(prec.getDriftTimeWindowUpperOffset() - 0.8) < 1e-5

    # Test charge
    prec.setCharge(2)
    assert prec.getCharge() == 2

    # Test possible charge states
    charges = [2,3,4]
    prec.setPossibleChargeStates(charges)
    assert list(prec.getPossibleChargeStates()) == charges

    # Test uncharged mass calculation
    mz = 200.0
    prec.setMZ(mz)
    charge = 2
    prec.setCharge(charge)
    expected_mass = mz * charge - charge * 1.007276466879  # mass of proton
    assert abs(prec.getUnchargedMass() - expected_mass) < 1e-5


@report
def testProcessingAction():
    """
    @tests: ProcessingAction
     ProcessingAction.ALIGNMENT
     ProcessingAction.BASELINE_REDUCTION
     ProcessingAction.CALIBRATION
     ProcessingAction.CHARGE_CALCULATION
     ProcessingAction.CHARGE_DECONVOLUTION
     ProcessingAction.CONVERSION_DTA
     ProcessingAction.CONVERSION_MZDATA
     ProcessingAction.CONVERSION_MZML
     ProcessingAction.CONVERSION_MZXML
     ProcessingAction.DATA_PROCESSING
     ProcessingAction.DEISOTOPING
     ProcessingAction.FEATURE_GROUPING
     ProcessingAction.FILTERING
     ProcessingAction.FORMAT_CONVERSION
     ProcessingAction.IDENTIFICATION
     ProcessingAction.IDENTIFICATION_MAPPING
     ProcessingAction.NORMALIZATION
     ProcessingAction.PEAK_PICKING
     ProcessingAction.PRECURSOR_RECALCULATION
     ProcessingAction.QUANTITATION
     ProcessingAction.SIZE_OF_PROCESSINGACTION
     ProcessingAction.SMOOTHING

    """
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.ALIGNMENT, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.BASELINE_REDUCTION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.CALIBRATION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.CHARGE_CALCULATION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.CHARGE_DECONVOLUTION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.CONVERSION_DTA, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.CONVERSION_MZDATA, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.CONVERSION_MZML, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.CONVERSION_MZXML, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.DATA_PROCESSING, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.DEISOTOPING, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.FEATURE_GROUPING, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.FILTERING, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.FORMAT_CONVERSION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.IDENTIFICATION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.IDENTIFICATION_MAPPING, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.NORMALIZATION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.PEAK_PICKING, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.PRECURSOR_RECALCULATION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.QUANTITATION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.SIZE_OF_PROCESSINGACTION, int)
    assert isinstance(pyopenms.DataProcessing.ProcessingAction.SMOOTHING, int)


@report
def testProduct():
    """
    @tests: Product
     Product.__init__
     Product.getIsolationWindowLowerOffset
     Product.getIsolationWindowUpperOffset
     Product.getMZ
     Product.setIsolationWindowLowerOffset
     Product.setIsolationWindowUpperOffset
     Product.setMZ
     Product.__eq__
     Product.__ge__
     Product.__gt__
     Product.__le__
     Product.__lt__
     Product.__ne__
    """
    p = pyopenms.Product()
    p.setMZ(12.0)
    p.setIsolationWindowLowerOffset(10.0)
    p.setIsolationWindowUpperOffset(15.0)
    assert p.getMZ() == 12.0
    assert p.getIsolationWindowLowerOffset() == 10.0
    assert p.getIsolationWindowUpperOffset() == 15.0

    assert p == p
    assert not p != p


@report
def testType():
    """
    @tests: Type
     Type.CONSENSUSXML
     Type.DTA
     Type.DTA2D
     Type.EDTA
     Type.FASTA
     Type.FEATUREXML
     Type.GELML
     Type.HARDKLOER
     Type.IDXML
     Type.INI
     Type.KROENIK
     Type.MASCOTXML
     Type.MGF
     Type.MS2
     Type.MSP
     Type.MZDATA
     Type.MZIDENTML
     Type.MZML
     Type.MZXML
     Type.OMSSAXML
     Type.PEPLIST
     Type.PEPXML
     Type.PNG
     Type.PROTXML
     Type.SIZE_OF_TYPE
     Type.TOPPAS
     Type.TRAML
     Type.TRANSFORMATIONXML
     Type.TSV
     Type.UNKNOWN
     Type.XMASS
    """
    for ti in  [
      pyopenms.FileType.CONSENSUSXML
     ,pyopenms.FileType.DTA
     ,pyopenms.FileType.DTA2D
     ,pyopenms.FileType.EDTA
     ,pyopenms.FileType.FASTA
     ,pyopenms.FileType.FEATUREXML
     ,pyopenms.FileType.GELML
     ,pyopenms.FileType.HARDKLOER
     ,pyopenms.FileType.IDXML
     ,pyopenms.FileType.INI
     ,pyopenms.FileType.KROENIK
     ,pyopenms.FileType.MASCOTXML
     ,pyopenms.FileType.MGF
     ,pyopenms.FileType.MS2
     ,pyopenms.FileType.MSP
     ,pyopenms.FileType.MZDATA
     ,pyopenms.FileType.MZIDENTML
     ,pyopenms.FileType.MZML
     ,pyopenms.FileType.MZXML
     ,pyopenms.FileType.OMSSAXML
     ,pyopenms.FileType.PEPLIST
     ,pyopenms.FileType.PEPXML
     ,pyopenms.FileType.PNG
     ,pyopenms.FileType.PROTXML
     ,pyopenms.FileType.SIZE_OF_TYPE
     ,pyopenms.FileType.TOPPAS
     ,pyopenms.FileType.TRAML
     ,pyopenms.FileType.TRANSFORMATIONXML
     ,pyopenms.FileType.TSV
     ,pyopenms.FileType.UNKNOWN
     ,pyopenms.FileType.XMASS]:
        assert isinstance(ti, int)


@report
def testVersion():
    """
    @tests: VersionDetails
     VersionDetails.__init__
     VersionDetails.create
     VersionDetails.version_major
     VersionDetails.version_minor
     VersionDetails.version_patch
     VersionDetails.__eq__
     VersionDetails.__ge__
     VersionDetails.__gt__
     VersionDetails.__le__
     VersionDetails.__lt__
     VersionDetails.__ne__
     VersionInfo.getRevision
     VersionInfo.getTime
     VersionInfo.getVersion
     __version__
    """
    _testStrOutput(pyopenms.VersionInfo.getVersion())
    _testStrOutput(pyopenms.VersionInfo.getRevision())
    _testStrOutput(pyopenms.VersionInfo.getTime())

    vd = pyopenms.VersionDetails.create("19.2.1")
    assert vd.version_major == 19
    assert vd.version_minor == 2
    assert vd.version_patch == 1

    vd = pyopenms.VersionDetails.create("19.2.1-alpha")
    assert vd.version_major == 19
    assert vd.version_minor == 2
    assert vd.version_patch == 1
    assert vd.pre_release_identifier == "alpha"

    assert vd == vd
    assert not vd < vd
    assert not vd > vd

    assert isinstance(pyopenms.__version__, str)


@report
def testAttachment():
    """
    @tests: Attachment
     Attachment.__init__
    """
    inst = pyopenms.Attachment()

    assert inst.name is not None
    assert inst.value is not None
    assert inst.cvRef is not None
    assert inst.cvAcc is not None
    assert inst.unitRef is not None
    assert inst.unitAcc is not None
    assert inst.binary is not None
    assert inst.qualityRef is not None
    assert inst.colTypes is not None
    assert inst.tableRows  is not None

    assert inst.toXMLString is not None
    assert inst.toCSVString is not None

    inst.name = "test"
    inst.value = "test"
    inst.cvRef = "test"
    inst.cvAcc = "test"
    inst.unitRef = "test"
    inst.unitAcc = "test"
    inst.binary = "test"
    inst.qualityRef = "test"
    inst.colTypes = [ b"test", b"test2"]
    inst.tableRows = [ [b"test", b"test2"], [b"otherTest"] ]

    assert inst.tableRows[1][0] == b"otherTest"


@report
def testRNaseDB():
    """
    @tests: RNaseDB
        const DigestionEnzymeRNA* getEnzyme(const String& name) except + nogil 
        const DigestionEnzymeRNA* getEnzymeByRegEx(const String& cleavage_regex) except + nogil 
        void getAllNames(libcpp_vector[ String ]& all_names) except + nogil 
        bool hasEnzyme(const String& name) except + nogil 
        bool hasRegEx(const String& cleavage_regex) except + nogil 
     """
    db = pyopenms.RNaseDB()
    names = []
    db.getAllNames(names)

    e = db.getEnzyme("RNase_T1")
    assert e.getThreePrimeGain() == u'p'

    assert db.hasEnzyme("RNase_T1")
    

@report
def testRibonucleotideDB():
    """
    @tests: RibonucleotideDB
    """
    r = pyopenms.RibonucleotideDB()

    uridine = r.getRibonucleotide(b"U")

    assert uridine.getName() == u'uridine'
    assert uridine.getCode() == u'U'
    assert uridine.getFormula().toString() == u'C9H12N2O6'
    assert uridine.isModified() == False


@report
def testRibonucleotide():
    """
    @tests: Ribonucleotide
    """
    r = pyopenms.Ribonucleotide()

    assert not r.isModified()

    r.setHTMLCode("test")
    assert r.getHTMLCode() == "test"

    r.setOrigin(b"A")
    assert r.getOrigin() == "A"

    r.setNewCode(b"A")
    assert r.getNewCode() == "A"


@report
def testRNaseDigestion():
    """
    @tests: RNaseDigestion
     """

    dig = pyopenms.RNaseDigestion()
    dig.setEnzyme("RNase_T1")
    assert dig.getEnzymeName() == "RNase_T1"

    oligo = pyopenms.NASequence.fromString("pAUGUCGCAG");

    result = []
    dig.digest(oligo, result)
    assert len(result) == 3


@report
def testNASequence():
    """
    @tests: NASequence
     """

    oligo = pyopenms.NASequence.fromString("pAUGUCGCAG");

    assert oligo.size() == 9
    seq_formula = oligo.getFormula()
    seq_formula.toString() == u'C86H108N35O64P9'

    oligo_mod = pyopenms.NASequence.fromString("A[m1A][Gm]A")
    seq_formula = oligo_mod.getFormula()
    seq_formula.toString() == u'C42H53N20O23P3'

    for r in oligo:
        pass

    assert oligo_mod[1].isModified()

    charge = 2
    oligo_mod.getMonoWeight(pyopenms.NASequence.NASFragmentType.WIon, charge)
    oligo_mod.getFormula(pyopenms.NASequence.NASFragmentType.WIon, charge)

    # Test __repr__ and __str__ methods
    # Test __repr__ method
    na_repr = pyopenms.NASequence.fromString("ACGU")
    repr_str = repr(na_repr)
    assert "NASequence(" in repr_str
    assert "sequence=" in repr_str
    assert "length=" in repr_str
    assert "mono_mass=" in repr_str
    
    # Test __str__ method returns the sequence string
    str_str = str(na_repr)
    assert str_str == "ACGU"

    # Test __len__ method
    assert len(na_repr) == 4


@report
def testExperimentalDesign():
    """
    @tests: ExperimentalDesign
     ExperimentalDesign.__init__
     ExperimentalDesign.getNumberOfSamples() == 8
     ExperimentalDesign.getNumberOfFractions() == 3
     ExperimentalDesign.getNumberOfLabels() == 4
     ExperimentalDesign.getNumberOfMSFiles() == 6
     ExperimentalDesign.getNumberOfFractionGroups() == 2
     ExperimentalDesign.getSample(1, 1) == 0
     ExperimentalDesign.getSample(2, 4) == 7
     ExperimentalDesign.isFractionated()
     ExperimentalDesign.sameNrOfMSFilesPerFraction()

     ExperimentalDesignFile.__init__
     ExperimentalDesignFile.load
     """
    f = pyopenms.ExperimentalDesignFile()
    fourplex_fractionated_design = pyopenms.ExperimentalDesign()
    ed_dirname = os.path.dirname(os.path.abspath(__file__))
    ed_filename = os.path.join(ed_dirname, "ExperimentalDesign_input_2.tsv").encode()
    fourplex_fractionated_design = pyopenms.ExperimentalDesignFile.load(ed_filename, False)
    assert fourplex_fractionated_design.getNumberOfSamples() == 8
    assert fourplex_fractionated_design.getNumberOfFractions() == 3
    assert fourplex_fractionated_design.getNumberOfLabels() == 4
    assert fourplex_fractionated_design.getNumberOfMSFiles() == 6
    assert fourplex_fractionated_design.getNumberOfFractionGroups() == 2
    assert fourplex_fractionated_design.getSample(1, 1) == 0
    assert fourplex_fractionated_design.getSample(2, 4) == 7
    assert fourplex_fractionated_design.isFractionated()
    assert fourplex_fractionated_design.sameNrOfMSFilesPerFraction()
 

@report
def testString():
    """
    @tests: String
     String.__init__
     String.toString
     String.c_str
    """
    pystr = pyopenms.String()
    pystr = pyopenms.String("blah")
    assert (pystr.toString() == "blah")
    pystr = pyopenms.String("blah")
    assert (pystr.toString() == "blah")
    pystr = pyopenms.String(u"blah")
    assert (pystr.toString() == "blah")
    pystr = pyopenms.String(pystr)
    assert (pystr.toString() == "blah")
    assert (len(pystr.toString()) == 4)
    cstr = pystr.c_str()

    # Printing should work ...
    print(cstr)
    print(pystr)
    print(pystr.toString())
    assert (pystr.toString() == "blah")

    pystr = pyopenms.String("bläh")
    assert (pystr.toString() == u"bläh")
    pystr = pyopenms.String("bläh")
    pystr = pyopenms.String(u"bläh")
    assert (pystr.toString() == u"bläh")
    pystr = pyopenms.String(pystr)
    assert (pystr.toString() == u"bläh")
    cstr = pystr.c_str()


if __name__ == "__main__":
    # Run all tests in this module
    import sys
    module = sys.modules[__name__]
    
    test_functions = [getattr(module, name) for name in dir(module) 
                     if name.startswith('test') and callable(getattr(module, name))]
    
    print(f"Running {len(test_functions)} miscellaneous tests...")
    for test_func in test_functions:
        try:
            test_func()
            print(f"✅ {test_func.__name__}")
        except Exception as e:
            print(f"❌ {test_func.__name__}: {e}")
