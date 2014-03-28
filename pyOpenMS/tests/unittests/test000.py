import pdb
import pyopenms

print "IMPORTED ", pyopenms.__file__


from functools import wraps


def report(f):
    @wraps(f)
    def wrapper(*a, **kw):
        print "run ", f.__name__
        f(*a, **kw)
    return wrapper



@report
def _testMetaInfoInterface(what):

    #void getKeys(libcpp_vector[String] & keys)
    #void getKeys(libcpp_vector[unsigned int] & keys)
    #DataValue getMetaValue(unsigned int) nogil except +
    #DataValue getMetaValue(String) nogil except +
    #void setMetaValue(unsigned int, DataValue) nogil except +
    #void setMetaValue(String, DataValue) nogil except +
    #bool metaValueExists(String) nogil except +
    #bool metaValueExists(unsigned int) nogil except +
    #void removeMetaValue(String) nogil except +
    #void removeMetaValue(unsigned int) nogil except +

    what.setMetaValue("key", 42)
    what.setMetaValue("key2", 42)

    keys = []
    what.getKeys(keys)
    keys = [0]
    what.getKeys(keys)
    assert len(keys) and all(isinstance(k, (long, int)) for k in keys)
    assert what.getMetaValue(keys[0]) == 42
    keys = [""]
    what.getKeys(keys)
    assert len(keys) and all(isinstance(k, str) for k in keys)

    assert what.getMetaValue(keys[0]) == 42

    assert what.metaValueExists("key")
    what.removeMetaValue("key")

    what.setMetaValue(1024, 42)

    keys = []
    what.getKeys(keys)
    keys = [0]
    what.getKeys(keys)
    assert len(keys) and all(isinstance(k, (long, int)) for k in keys)
    assert what.getMetaValue(keys[0]) == 42
    keys = [""]
    what.getKeys(keys)
    assert len(keys) and all(isinstance(k, str) for k in keys)

    assert what.getMetaValue(keys[0]) == 42

    what.setMetaValue("key", 42)
    what.setMetaValue("key2", 42)

    assert what.metaValueExists("key")
    what.removeMetaValue("key")
    keys = []
    what.getKeys(keys)
    assert len(keys) == 1
    what.removeMetaValue("key2")
    keys = []
    what.getKeys(keys)
    assert len(keys) == 0


    what.clearMetaInfo()
    keys = []
    what.getKeys(keys)
    assert len(keys) == 0


@report
def _testUniqueIdInterface(what):

    assert what.hasInvalidUniqueId()
    assert not what.hasValidUniqueId()
    assert what.ensureUniqueId()
    assert isinstance(what.getUniqueId(), (int, long))
    assert what.getUniqueId() > 0
    assert not what.hasInvalidUniqueId()
    assert what.hasValidUniqueId()

    what.clearUniqueId()
    assert what.getUniqueId() == 0
    assert what.hasInvalidUniqueId()
    assert not what.hasValidUniqueId()

    assert what.ensureUniqueId()
    assert isinstance(what.getUniqueId(), (int, long))
    assert what.getUniqueId() > 0
    assert not what.hasInvalidUniqueId()
    assert what.hasValidUniqueId()

    what.setUniqueId(1234)
    assert  what.getUniqueId() == 1234


def _testProgressLogger(ff):
    """
    @tests:
     ProgressLogger.__init__
     ProgressLogger.endProgress
     ProgressLogger.getLogType
     ProgressLogger.setLogType
     ProgressLogger.setProgress
     ProgressLogger.startProgress
    """

    ff.setLogType(pyopenms.LogType.NONE)
    assert ff.getLogType() == pyopenms.LogType.NONE
    ff.startProgress(0, 3, "label")
    ff.setProgress(0)
    ff.setProgress(1)
    ff.setProgress(2)
    ff.setProgress(3)
    ff.endProgress()

@report
def testSpectrumAlignment():
    """
    @tests:
        SpectrumAlignment.__init__
        SpectrumAlignment.getSpectrumAlignment
    """
    # test existence of some methods
    pyopenms.SpectrumAlignment
    pyopenms.SpectrumAlignment.__init__
    pyopenms.SpectrumAlignment.getDefaults
    pyopenms.SpectrumAlignment.getParameters
    pyopenms.SpectrumAlignment.setParameters

    spec = pyopenms.MSSpectrum()
    p = pyopenms.Peak1D()
    p.setMZ(1000.0)
    p.setIntensity(200.0)
    spec.push_back(p)
    p.setMZ(2000.0)
    p.setIntensity(200.0)
    spec.push_back(p)

    rich_spec = pyopenms.RichMSSpectrum()
    p = pyopenms.RichPeak1D()
    p.setMZ(1000.001)
    p.setIntensity(200.0)
    rich_spec.push_back(p)
    p.setMZ(2000.001)
    p.setIntensity(200.0)
    rich_spec.push_back(p)
    p.setMZ(3000.001)
    p.setIntensity(200.0)
    rich_spec.push_back(p)

    aligner = pyopenms.SpectrumAlignment()
    result = []

    aligner.getSpectrumAlignment(result, spec, spec)
    assert result == [ (0,0), (1,1) ], result
    aligner.getSpectrumAlignment(result, rich_spec, spec)
    assert result == [ (0,0), (1,1) ], result
    aligner.getSpectrumAlignment(result, spec, rich_spec)
    assert result == [ (0,0), (1,1) ], result
    aligner.getSpectrumAlignment(result, rich_spec, rich_spec)
    assert result == [ (0,0), (1,1), (2,2) ], result

    aligner = pyopenms.SpectrumAlignmentScore()
    assert isinstance(aligner(spec), float)
    assert isinstance(aligner(rich_spec), float)

    assert isinstance(aligner(spec, rich_spec), float)
    assert isinstance(aligner(rich_spec, spec), float)

    assert isinstance(aligner(spec, spec), float)
    assert isinstance(aligner(rich_spec, rich_spec), float)


@report
def testAASequence():
    """
    @tests:
     AASequence.__init__
     AASequence.__add__
     AASequence.__radd__
     AASequence.__iadd__
     AASequence.getCTerminalModification
     AASequence.getNTerminalModification
     AASequence.setCTerminalModification
     AASequence.setModification
     AASequence.setNTerminalModification
     AASequence.setStringSequence
     AASequence.toString
     AASequence.toUnmodifiedString
    """
    aas = pyopenms.AASequence()

    aas + aas
    aas += aas

    aas.__doc__
    aas = pyopenms.AASequence("DFPIANGER")
    assert aas.getCTerminalModification() == ""
    assert aas.getNTerminalModification() == ""
    aas.setCTerminalModification("")
    aas.setNTerminalModification("")
    aas.setStringSequence("")
    assert aas.toString() == ""
    assert aas.toUnmodifiedString() == ""

@report
def testElement():
    """
    @tests:
     Element.__init__
     Element.setAtomicNumber
     Element.getAtomicNumber
     Element.setAverageWeight
     Element.getAverageWeight
     Element.setMonoWeight
     Element.getMonoWeight
     Element.setIsotopeDistribution
     Element.getIsotopeDistribution
     Element.setName
     Element.getName
     Element.setSymbol
     Element.getSymbol
    """
    ins = pyopenms.Element()

    ins.setAtomicNumber(6)
    ins.getAtomicNumber()
    ins.setAverageWeight(12.011)
    ins.getAverageWeight()
    ins.setMonoWeight(12)
    ins.getMonoWeight()
    iso = pyopenms.IsotopeDistribution()
    ins.setIsotopeDistribution(iso)
    ins.getIsotopeDistribution()
    ins.setName("Carbon")
    ins.getName()
    ins.setSymbol("C")
    ins.getSymbol()

@report
def testResidue():
    """
    @tests:
     Residue.__init__
    """
    ins = pyopenms.Residue()

    pyopenms.Residue.ResidueType.Full
    pyopenms.Residue.ResidueType.Internal
    pyopenms.Residue.ResidueType.NTerminal
    pyopenms.Residue.ResidueType.CTerminal
    pyopenms.Residue.ResidueType.AIon
    pyopenms.Residue.ResidueType.BIon
    pyopenms.Residue.ResidueType.CIonMinusOne
    pyopenms.Residue.ResidueType.CIon
    pyopenms.Residue.ResidueType.CIonPlusOne
    pyopenms.Residue.ResidueType.CIonPlusTwo
    pyopenms.Residue.ResidueType.XIon
    pyopenms.Residue.ResidueType.YIon
    pyopenms.Residue.ResidueType.ZIonMinusOne
    pyopenms.Residue.ResidueType.ZIon
    pyopenms.Residue.ResidueType.ZIonPlusOne
    pyopenms.Residue.ResidueType.ZIonPlusTwo
    pyopenms.Residue.ResidueType.SizeOfResidueType

@report
def testIsotopeDistribution():
    """
    @tests:
     IsotopeDistribution.__init__
    """
    ins = pyopenms.IsotopeDistribution()

    ins.setMaxIsotope(5)
    ins.getMaxIsotope()
    ins.getMax()
    ins.getMin()
    ins.size()
    ins.clear()
    ins.estimateFromPeptideWeight(500)
    ins.renormalize()
    ins.trimLeft(6.0)
    ins.trimRight(8.0)

@report
def testEmpiricalFormula():
    """
    @tests:
     EmpiricalFormula.__init__
     EmpiricalFormula.getMonoWeight
     EmpiricalFormula.getAverageWeight
     EmpiricalFormula.getIsotopeDistribution
     EmpiricalFormula.getNumberOfAtoms
     EmpiricalFormula.setCharge
     EmpiricalFormula.getCharge
     EmpiricalFormula.toString
     EmpiricalFormula.isEmpty
     EmpiricalFormula.isCharged
     EmpiricalFormula.hasElement
     EmpiricalFormula.hasElement
    """
    ins = pyopenms.EmpiricalFormula()

    ins.getMonoWeight()
    ins.getAverageWeight()
    ins.getIsotopeDistribution(1)
    # ins.getNumberOf(0)
    # ins.getNumberOf("test")
    ins.getNumberOfAtoms()
    ins.setCharge(2)
    ins.getCharge()
    ins.toString()
    ins.isEmpty()
    ins.isCharged()
    ins.hasElement("C")
    ins.hasElement(6)

@report
def testIdentificationHit():
    """
    @tests:
     IdentificationHit.__init__
    """
    f = pyopenms.IdentificationHit()
    _testMetaInfoInterface(f)

    assert pyopenms.IdentificationHit().setId is not None
    assert pyopenms.IdentificationHit().getId is not None
    assert pyopenms.IdentificationHit().setCharge is not None
    assert pyopenms.IdentificationHit().getCharge is not None
    assert pyopenms.IdentificationHit().setCalculatedMassToCharge is not None
    assert pyopenms.IdentificationHit().getCalculatedMassToCharge is not None
    assert pyopenms.IdentificationHit().setExperimentalMassToCharge is not None
    assert pyopenms.IdentificationHit().getExperimentalMassToCharge is not None
    assert pyopenms.IdentificationHit().setName is not None
    assert pyopenms.IdentificationHit().getName is not None
    assert pyopenms.IdentificationHit().setPassThreshold is not None
    assert pyopenms.IdentificationHit().getPassThreshold is not None
    assert pyopenms.IdentificationHit().setRank is not None
    assert pyopenms.IdentificationHit().getRank is not None

    f.setId("test_id")
    assert f.getId() == "test_id"

    f.setCharge(5)
    assert f.getCharge() == 5

    f.setCalculatedMassToCharge(5.0)
    assert f.getCalculatedMassToCharge() == 5.0

    f.setExperimentalMassToCharge(5.0)
    assert f.getExperimentalMassToCharge() == 5.0

    f.setName("test")
    assert f.getName() == "test"

    f.setPassThreshold(True)
    assert f.getPassThreshold() == True

    f.setRank(42)
    assert f.getRank() == 42

@report
def testSpectrumIdentification():
    """
    @tests:
     SpectrumIdentification.__init__
    """
    f = pyopenms.SpectrumIdentification()
    _testMetaInfoInterface(f)

    assert pyopenms.SpectrumIdentification().setHits is not None
    assert pyopenms.SpectrumIdentification().addHit is not None
    assert pyopenms.SpectrumIdentification().getHits is not None

    hit = pyopenms.IdentificationHit()
    hit.setName("test1")
    f.addHit(hit)
    hit = pyopenms.IdentificationHit()
    hit.setName("test2")
    f.addHit(hit)
    all_hits = f.getHits()
    assert len(all_hits) == 2
    assert "test1" in [h.getName() for h in all_hits]
    assert "test2" in [h.getName() for h in all_hits]

@report
def testIdentification():
    """
    @tests:
     Identification.__init__
    """
    f = pyopenms.Identification()
    _testMetaInfoInterface(f)

    assert pyopenms.Identification().setCreationDate is not None
    assert pyopenms.Identification().getCreationDate is not None
    assert pyopenms.Identification().setSpectrumIdentifications is not None
    assert pyopenms.Identification().addSpectrumIdentification is not None
    assert pyopenms.Identification().getSpectrumIdentifications is not None

    id1 = pyopenms.SpectrumIdentification()
    f.addSpectrumIdentification(id1)
    assert len(f.getSpectrumIdentifications()) == 1
    id2 = pyopenms.SpectrumIdentification()
    f.addSpectrumIdentification(id2)
    assert len(f.getSpectrumIdentifications()) == 2

@report
def testModificationDefinitionsSet():
    """
    @tests:
     ModificationDefinitionsSet.__init__
    """
    empty = pyopenms.ModificationDefinitionsSet()
    fixed = ["Carbamidomethyl"]
    variable = ["Oxidation"]
    full = pyopenms.ModificationDefinitionsSet(fixed, variable)

@report
def test_AcquisitionInfo():
    """
    @tests:
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
def test_BaseFeature():
    """
    @tests:
     BaseFeature.__init__
     BaseFeature.clearUniqueId
     BaseFeature.ensureUniqueId
     BaseFeature.getCharge
     BaseFeature.getKeys
     BaseFeature.getMetaValue
     BaseFeature.getQuality
     BaseFeature.getUniqueId
     BaseFeature.getWidth
     BaseFeature.hasInvalidUniqueId
     BaseFeature.hasValidUniqueId
     BaseFeature.metaValueExists
     BaseFeature.removeMetaValue
     BaseFeature.setCharge
     BaseFeature.setMetaValue
     BaseFeature.setQuality
     BaseFeature.setWidth
     BaseFeature.clearMetaInfo
     BaseFeature.setUniqueId
    """
    bf = pyopenms.BaseFeature()
    _testMetaInfoInterface(bf)
    _testUniqueIdInterface(bf)
    bf.clearUniqueId()
    assert bf.ensureUniqueId()
    assert bf.getCharge() == 0
    assert isinstance(bf.getQuality(), float)
    assert isinstance(bf.getUniqueId(), (long, int))
    assert isinstance(bf.getWidth(), float)

    assert not bf.hasInvalidUniqueId()
    assert bf.hasValidUniqueId()


    _testMetaInfoInterface(bf)
    bf.setCharge(1)
    bf.setQuality(0.0)
    bf.setWidth(1.0)

@report
def test_AnnotationState():
    """
    @tests:
     AnnotationState.__init__
    """
    state = pyopenms.AnnotationState()

    assert state.FEATURE_ID_NONE is not None
    assert state.FEATURE_ID_SINGLE is not None
    assert state.FEATURE_ID_MULTIPLE_SAME is not None
    assert state.FEATURE_ID_MULTIPLE_DIVERGENT is not None
    assert state.SIZE_OF_ANNOTATIONSTATE is not None

@report
def testChecksumType():
    """
    @tests:
     ChecksumType.MD5
     ChecksumType.SHA1
     ChecksumType.SIZE_OF_CHECKSUMTYPE
     ChecksumType.UNKNOWN_CHECKSUM
    """
    assert isinstance(pyopenms.ChecksumType.MD5, int)
    assert isinstance(pyopenms.ChecksumType.SHA1, int)
    assert isinstance(pyopenms.ChecksumType.SIZE_OF_CHECKSUMTYPE, int)
    assert isinstance(pyopenms.ChecksumType.UNKNOWN_CHECKSUM, int)


@report
def testChromatogramPeak():
    """
    @tests:
     ChromatogramPeak.__init__
     ChromatogramPeak.__eq__
     ChromatogramPeak.__ge__
     ChromatogramPeak.__gt__
     ChromatogramPeak.__le__
     ChromatogramPeak.__lt__
     ChromatogramPeak.__ne__
     ChromatogramPeak.getIntensity
     ChromatogramPeak.getRT
     ChromatogramPeak.setIntensity
     ChromatogramPeak.setRT
    """
    p = pyopenms.ChromatogramPeak()
    assert p == p
    assert not p != p


    p.setIntensity(12.0)
    p.setRT(34.0)
    assert p.getIntensity() == 12.0
    assert p.getRT() == 34.0



@report
def testChromatogramToosl():
    """
    @tests:
     ChromatogramTools.__init__
     ChromatogramTools.convertChromatogramsToSpectra
     ChromatogramTools.convertSpectraToChromatograms
    """
    pyopenms.ChromatogramTools()
    pyopenms.ChromatogramTools.convertChromatogramsToSpectra
    pyopenms.ChromatogramTools.convertSpectraToChromatograms


@report
def testConsensusFeature():
    """
    @tests:
     ConsensusFeature.__eq__
     ConsensusFeature.__ge__
     ConsensusFeature.__gt__
     ConsensusFeature.__le__
     ConsensusFeature.__lt__
     ConsensusFeature.__ne__
     ConsensusFeature.__init__
     ConsensusFeature.clearUniqueId
     ConsensusFeature.computeConsensus
     ConsensusFeature.computeDechargeConsensus
     ConsensusFeature.computeMonoisotopicConsensus
     ConsensusFeature.ensureUniqueId
     ConsensusFeature.getCharge
     ConsensusFeature.getKeys
     ConsensusFeature.getMetaValue
     ConsensusFeature.getQuality
     ConsensusFeature.getUniqueId
     ConsensusFeature.getWidth
     ConsensusFeature.hasInvalidUniqueId
     ConsensusFeature.hasValidUniqueId
     ConsensusFeature.insert
     ConsensusFeature.metaValueExists
     ConsensusFeature.removeMetaValue
     ConsensusFeature.setCharge
     ConsensusFeature.setMetaValue
     ConsensusFeature.setQuality
     ConsensusFeature.setWidth
     ConsensusFeature.clearMetaInfo
     ConsensusFeature.setUniqueId
     ConsensusFeature.size

     ConsensusFeature.getPeptideIdentifications
     ConsensusFeature.setPeptideIdentifications
    """


    f = pyopenms.ConsensusFeature()
    _testUniqueIdInterface(f)
    _testMetaInfoInterface(f)

    f.setCharge(1)
    f.setQuality(2.0)
    f.setWidth(4.0)
    assert f.getCharge() == 1
    assert f.getQuality() == 2.0
    assert f.getWidth() == 4.0

    f.insert(0, pyopenms.Peak2D(), 1)
    f.insert(1, pyopenms.BaseFeature())
    f.insert(2, pyopenms.ConsensusFeature())

    f.computeConsensus()
    f.computeDechargeConsensus
    f.computeMonoisotopicConsensus()

    assert f.size() >= 0

    p = f.getPeptideIdentifications()
    f.setPeptideIdentifications(p)


@report
def testConsensusMap():
    """
    @tests:
     ConsensusMap.__eq__
     ConsensusMap.__ge__
     ConsensusMap.__gt__
     ConsensusMap.__init__
     ConsensusMap.__iter__
     ConsensusMap.__le__
     ConsensusMap.__lt__
     ConsensusMap.__ne__
     ConsensusMap.clear
     ConsensusMap.clearUniqueId
     ConsensusMap.ensureUniqueId
     ConsensusMap.getDataProcessing
     ConsensusMap.getFileDescriptions
     ConsensusMap.getProteinIdentifications
     ConsensusMap.getUnassignedPeptideIdentifications
     ConsensusMap.getUniqueId
     ConsensusMap.hasInvalidUniqueId
     ConsensusMap.hasValidUniqueId
     ConsensusMap.setDataProcessing
     ConsensusMap.setFileDescriptions
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

    m.clear()
    m.clearUniqueId()
    m.ensureUniqueId()
    m.getDataProcessing()
    m.getFileDescriptions()
    m.getProteinIdentifications()
    m.getUnassignedPeptideIdentifications()
    m.getUniqueId()
    m.hasInvalidUniqueId()
    m.hasValidUniqueId()
    m.setDataProcessing
    m.setFileDescriptions
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
    m.updateRanges()

    assert isinstance(m.getMin()[0], float)
    assert isinstance(m.getMin()[0], float)
    assert isinstance(m.getMax()[1], float)
    assert isinstance(m.getMax()[1], float)
    assert isinstance(m.getMinInt(), float)
    assert isinstance(m.getMaxInt(), float)

    m.getIdentifier()
    m.getLoadedFileType()
    m.getLoadedFilePath()

    assert m == m
    assert not m != m

@report
def testConsensusXMLFile():
    """
    @tests:
     ConsensusXMLFile.__init__
     ConsensusXMLFile.getOptions
     ConsensusXMLFile.load
     ConsensusXMLFile.store
    """
    f = pyopenms.ConsensusXMLFile()
    f.getOptions()
    assert f.load is not None
    assert f.store is not None

@report
def testXTandemXMLFile():
    """
    @tests:
     XTandemXMLFile.__init__
     XTandemXMLFile.load
     XTandemXMLFile.setModificationDefinitionsSet
    """
    f = pyopenms.XTandemXMLFile()
    assert f.load is not None
    assert f.setModificationDefinitionsSet is not None

@report
def testSignalToNoiseEstimatorMedian():
    """
    @tests:
     SignalToNoiseEstimatorMedian.__init__
    """
    f = pyopenms.SignalToNoiseEstimatorMedian()
    assert f.init is not None
    assert f.getSignalToNoise is not None

@report
def testSignalToNoiseEstimatorMedianChrom():
    """
    @tests:
     SignalToNoiseEstimatorMedianChrom.__init__
    """
    f = pyopenms.SignalToNoiseEstimatorMedianChrom()
    assert f.init is not None
    assert f.getSignalToNoise is not None

@report
def testConvexHull2D():
    """
    @tests:
     ConvexHull2D.__eq__
     ConvexHull2D.__ge__
     ConvexHull2D.__gt__
     ConvexHull2D.__init__
     ConvexHull2D.__le__
     ConvexHull2D.__lt__
     ConvexHull2D.__ne__
     ConvexHull2D.clear
     """
    ch = pyopenms.ConvexHull2D()
    ch.clear()
    assert ch == ch


@report
def testDataProcessing(dp=pyopenms.DataProcessing()):

    """
    @tests:
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

    _testMetaInfoInterface(dp)

    assert dp == dp
    assert not dp != dp

    assert isinstance(dp.getCompletionTime().getDate(), str)
    assert isinstance(dp.getCompletionTime().getTime(), str)
    dp.clearMetaInfo()
    k = []
    dp.getKeys(k)
    assert k == []
    dp.getMetaValue
    ac = dp.getProcessingActions()
    assert ac == set(())
    dp.setProcessingActions(ac)
    assert isinstance(dp.getSoftware().getName(), str)
    assert isinstance(dp.getSoftware().getVersion(), str)
    dp.isMetaEmpty()
    dp.metaValueExists
    dp.removeMetaValue
    dp.setCompletionTime(pyopenms.DateTime.now())
    s = dp.getSoftware()
    s.setName("pyopenms")
    dp.setSoftware(s)

    assert dp.getSoftware().getName() == "pyopenms"


@report
def testDataType():
    """
    @tests:
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
    @tests:
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

    a = pyopenms.DataValue(["1.0"])
    assert not a.isEmpty()
    assert a.toStringList() == ["1.0"]
    assert a.valueType() == pyopenms.DataType.STRING_LIST

@report
def testAdduct():
    """
    @tests:
     Adduct.__init__
    """
    a = pyopenms.Adduct()

@report
def testGaussFitter():
    """
    @tests:
     GaussFitter.__init__
    """
    ins = pyopenms.GaussFitter()

@report
def testGaussFitResult():
    """
    @tests:
     GaussFitResult.__init__
    """
    ins = pyopenms.GaussFitResult(0.0, 0.0, 0.0)
    ins.A = 5.0
    ins.x0 = 5.0
    ins.sigma = 5.0

@report
def testChargePair():
    """
    @tests:
     ChargePair.__init__
    """
    a = pyopenms.ChargePair()

@report
def testCompomer():
    """
    @tests:
     Compomer.__init__
    """
    a = pyopenms.Compomer()

@report
def testCVMappings():
    """
    @tests:
     CVMappings.__init__
    """
    val = pyopenms.CVMappings()

@report
def testCVMappingFile():
    """
    @tests:
     CVMappingFile.__init__
    """
    val = pyopenms.CVMappingFile()

    assert pyopenms.CVMappingFile().load

@report
def testControlledVocabulary():
    """
    @tests:
     ControlledVocabulary.__init__
    """
    val = pyopenms.ControlledVocabulary()

    assert pyopenms.ControlledVocabulary().loadFromOBO

@report
def testSemanticValidator():
    """
    @tests:
     SemanticValidator.__init__
    """
    m = pyopenms.CVMappings()
    cv = pyopenms.ControlledVocabulary()

    val = pyopenms.SemanticValidator(m, cv)

    assert val.validate is not None
    assert val.setCheckTermValueTypes is not None
    assert val.setCheckUnits is not None


@report
def testDateTime():
    """
    @tests:
     DateTime.__init__
     DateTime.getDate
     DateTime.getTime
     DateTime.now
    """
    d = pyopenms.DateTime()
    assert isinstance( d.getDate(), str)
    assert isinstance( d.getTime(), str)
    d = pyopenms.DateTime.now()
    assert isinstance( d.getDate(), str)
    assert isinstance( d.getTime(), str)

    d.clear()
    d.set("01.01.2001 11:11:11")
    assert d.get() == "2001-01-01 11:11:11"

@report
def testFeature():
    """
    @tests:
     Feature.__init__
     Feature.clearUniqueId
     Feature.ensureUniqueId
     Feature.getCharge
     Feature.getIntensity
     Feature.getKeys
     Feature.getMZ
     Feature.getMetaValue
     Feature.getQuality
     Feature.getRT
     Feature.getUniqueId
     Feature.getWidth
     Feature.hasInvalidUniqueId
     Feature.hasValidUniqueId
     Feature.metaValueExists
     Feature.removeMetaValue
     Feature.setCharge
     Feature.setIntensity
     Feature.setMZ
     Feature.setMetaValue
     Feature.setQuality
     Feature.setRT
     Feature.setWidth
     Feature.__eq__
     Feature.__ge__
     Feature.__gt__
     Feature.__le__
     Feature.__lt__
     Feature.__ne__
     Feature.clearMetaInfo

     Feature.getConvexHulls
     Feature.getSubordinates
     Feature.setConvexHulls
     Feature.setSubordinates
     Feature.setUniqueId

     Feature.getPeptideIdentifications
     Feature.setPeptideIdentifications
    """
    f = pyopenms.Feature()
    _testMetaInfoInterface(f)
    _testUniqueIdInterface(f)

    f.setConvexHulls(f.getConvexHulls())
    f.setSubordinates(f.getSubordinates())
    f.setUniqueId(12345)

    assert f == f
    assert not f != f

    f.setCharge(-1)
    assert f.getCharge() == -1
    f.setIntensity(10.0)
    assert f.getIntensity() == 10.0
    f.setQuality(0, 20.0)
    assert f.getQuality(0) == 20.0
    f.setRT(30.0)
    assert f.getRT() == 30.0
    f.setWidth(40.0)
    assert f.getWidth() == 40.0

    p = f.getPeptideIdentifications()
    f.setPeptideIdentifications(p)


@report
def testFeatureFinder():
    """
    @tests:
     FeatureFinder.__init__
     FeatureFinder.endProgress
     FeatureFinder.getLogType
     FeatureFinder.getParameters
     FeatureFinder.run
     FeatureFinder.setLogType
     FeatureFinder.setProgress
     FeatureFinder.startProgress
    """
    ff = pyopenms.FeatureFinder()
    name = pyopenms.FeatureFinderAlgorithmPicked.getProductName()
    ff.run(name, pyopenms.MSExperiment(), pyopenms.FeatureMap() ,
            pyopenms.Param(), pyopenms.FeatureMap())

    _testProgressLogger(ff)

    p = ff.getParameters(name)
    _testParam(p)

@report
def testFeatureFileOptions():
    """
    @tests:
     FeatureFileOptions.__init__
     FeatureFileOptions.getLoadConvexHull
     FeatureFileOptions.getLoadSubordinates
     FeatureFileOptions.getMetadataOnly
     FeatureFileOptions.getSizeOnly
     FeatureFileOptions.setLoadConvexHull
     FeatureFileOptions.setLoadSubordinates
     FeatureFileOptions.setMetadataOnly
     FeatureFileOptions.setSizeOnly
    """

    fo = pyopenms.FeatureFileOptions()
    fo.getLoadConvexHull()
    fo.getLoadSubordinates()
    fo.getSizeOnly()
    assert fo.setLoadConvexHull is not None
    assert fo.setLoadSubordinates is not None
    assert fo.setMetadataOnly is not None
    assert fo.setSizeOnly is not None

@report
def _testParam(p):
    """
    @tests:
     Param.__init__
     Param.addTag
     Param.addTags
     Param.asDict
     Param.clearTags
     Param.copy
     Param.exists
     Param.get
     Param.getDescription
     Param.getEntry
     Param.getSectionDescription
     Param.getTags
     Param.getValue
     Param.hasTag
     Param.insert
     Param.setMaxFloat
     Param.setMaxInt
     Param.setMinFloat
     Param.setMinInt
     Param.setSectionDescription
     Param.setValidStrings
     Param.setValue
     Param.size
     Param.update
     Param.__eq__
     Param.__ge__
     Param.__gt__
     Param.__le__
     Param.__lt__
     Param.__ne__
     ParamEntry.__init__
     ParamEntry.description
     ParamEntry.isValid
     ParamEntry.max_float
     ParamEntry.max_int
     ParamEntry.min_float
     ParamEntry.min_int
     ParamEntry.name
     ParamEntry.tags
     ParamEntry.valid_strings
     ParamEntry.value
     ParamEntry.__eq__
     ParamEntry.__ge__
     ParamEntry.__gt__
     ParamEntry.__le__
     ParamEntry.__lt__
     ParamEntry.__ne__
    """

    assert p == p

    dd = p.asDict()
    assert len(dd) == p.size()
    assert isinstance(dd, dict)

    for k in p.keys():
        #value = p.getValue(k)
        value = p[k]
        p[k] = value
        p.update(p)
        p.update(p.asDict())
        assert p[k] == value
        desc  = p.getDescription(k)
        tags  = p.getTags(k)
        p.setValue(k, value, desc, tags)
        p.setValue(k, value, desc)
        assert p.exists(k)
        # only set the section description if there are actully two or more sections
        if len(k.split(":")) < 2: continue
        f = k.split(":")[0]
        p.setSectionDescription(f, k)
        assert p.getSectionDescription(f) == k

        assert p.get(k) is not None

    assert sorted(p.items()) == sorted((k, p[k]) for k in p.keys())
    assert sorted(p.values()) == sorted(p[k] for k in p.keys())


    assert not p.exists("asdflkj01231321321v")
    p.addTag(k, "a")
    p.addTags(k, ["b", "c"])
    assert sorted(p.getTags(k)) == ["a", "b", "c"]
    p.clearTags(k)
    assert p.getTags(k) == []

    pn = pyopenms.Param()
    pn.insert("master:", p)
    assert pn.exists("master:"+k)

    p1 = pn.copy("master:", True)
    assert p1 == p

    p1.update(p)
    p1.update(p,0)
    p1.update(p,1)
    p1.update(dd)

    p.setValidStrings
    p.setMinFloat
    p.setMaxFloat
    p.setMinInt
    p.setMaxInt
    ph = pyopenms.ParamXMLFile()
    ph.store("test.ini", p)
    p1 = pyopenms.Param()
    ph.load("test.ini", p1)
    assert p == p1

    e1 = p1.getEntry(k)
    for f in ["name", "description", "value", "tags", "valid_strings",
              "min_float", "max_float", "min_int", "max_int"]:
        assert getattr(e1, f) is not None

    assert e1 == e1

    assert p1.get("abcde", 7) == 7




@report
def testFeatureFinderAlgorithmPicked():
    """
    @tests:
     FeatureFinderAlgorithmPicked.__init__
     FeatureFinderAlgorithmPicked.getDefaults
     FeatureFinderAlgorithmPicked.getName
     FeatureFinderAlgorithmPicked.getParameters
     FeatureFinderAlgorithmPicked.getProductName
     FeatureFinderAlgorithmPicked.setName
     FeatureFinderAlgorithmPicked.setParameters
    """
    ff = pyopenms.FeatureFinderAlgorithmPicked()
    p = ff.getDefaults()
    _testParam(p)

    _testParam(ff.getParameters())

    assert ff.getName() == "FeatureFinderAlgorithm"
    assert pyopenms.FeatureFinderAlgorithmPicked.getProductName() == "centroided"

    ff.setParameters(pyopenms.Param())

    ff.setName("test")
    assert ff.getName() == "test"

@report
def testFeatureFinderAlgorithmSH():
    """
    @tests:
     FeatureFinderAlgorithmSH.__init__
     FeatureFinderAlgorithmSH.getDefaults
     FeatureFinderAlgorithmSH.getName
     FeatureFinderAlgorithmSH.getParameters
     FeatureFinderAlgorithmSH.getProductName
     FeatureFinderAlgorithmSH.setName
     FeatureFinderAlgorithmSH.setParameters
    """
    ff = pyopenms.FeatureFinderAlgorithmSH()
    p = ff.getDefaults()
    _testParam(p)

    # _testParam(ff.getParameters())

    assert ff.getName() == "FeatureFinderAlgorithm"
    assert pyopenms.FeatureFinderAlgorithmSH.getProductName() == "superhirn"

    ff.setParameters(pyopenms.Param())

    ff.setName("test")
    assert ff.getName() == "test"

@report
def testFeatureFinderAlgorithmIsotopeWavelet():
    """
    @tests:
     FeatureFinderAlgorithmIsotopeWavelet.__init__
     FeatureFinderAlgorithmIsotopeWavelet.getDefaults
     FeatureFinderAlgorithmIsotopeWavelet.getName
     FeatureFinderAlgorithmIsotopeWavelet.getParameters
     FeatureFinderAlgorithmIsotopeWavelet.getProductName
     FeatureFinderAlgorithmIsotopeWavelet.setName
     FeatureFinderAlgorithmIsotopeWavelet.setParameters
    """
    ff = pyopenms.FeatureFinderAlgorithmIsotopeWavelet()
    p = ff.getDefaults()
    _testParam(p)

    # _testParam(ff.getParameters())

    assert ff.getName() == "FeatureFinderAlgorithm"
    assert pyopenms.FeatureFinderAlgorithmIsotopeWavelet.getProductName() == "isotope_wavelet"

    ff.setParameters(pyopenms.Param())

    ff.setName("test")
    assert ff.getName() == "test"

@report
def testCompNovoIdentification():
    """
    @tests:
     CompNovoIdentification.__init__
    """
    ff = pyopenms.CompNovoIdentification()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.CompNovoIdentification().getIdentification is not None
    assert pyopenms.CompNovoIdentification().getIdentifications is not None

@report
def testCompNovoIdentificationCID():
    """
    @tests:
     CompNovoIdentificationCID.__init__
    """
    ff = pyopenms.CompNovoIdentificationCID()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.CompNovoIdentificationCID().getIdentification is not None
    assert pyopenms.CompNovoIdentificationCID().getIdentifications is not None

@report
def testExperimentalSettings():
    """
    @tests:
     ExperimentalSettings.__init__
    """
    ff = pyopenms.ExperimentalSettings()

@report
def testFeatureDeconvolution():
    """
    @tests:
     FeatureDeconvolution.__init__
    """
    ff = pyopenms.FeatureDeconvolution()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.FeatureDeconvolution().compute is not None

@report
def testInternalCalibration():
    """
    @tests:
     InternalCalibration.__init__
    """
    ff = pyopenms.InternalCalibration()
    p = ff.getDefaults()
    _testParam(p)

    # TODO 
    # assert pyopenms.InternalCalibration().compute is not None

@report
def testItraqConstants():
    """
    @tests:
    """
    constants = pyopenms.ItraqConstants()

    assert pyopenms.ITRAQ_TYPES.FOURPLEX is not None
    assert pyopenms.ITRAQ_TYPES.EIGHTPLEX is not None
    assert pyopenms.ITRAQ_TYPES.TMT_SIXPLEX is not None

    assert constants.getIsotopeMatrixAsStringList is not None
    assert constants.updateIsotopeMatrixFromStringList is not None
    assert constants.translateIsotopeMatrix is not None

@report
def testItraqChannelExtractor():
    """
    @tests:
     ItraqChannelExtractor.__init__
    """
    extractor = pyopenms.ItraqChannelExtractor()
    p = extractor.getDefaults()
    _testParam(p)


    # Note that using TMT_SIXPLEX will not work here
    p = extractor.getDefaults()
    pyopenms.ItraqChannelExtractor(pyopenms.ITRAQ_TYPES.FOURPLEX, p)
    assert pyopenms.ItraqChannelExtractor().run is not None

    assert extractor.getIsotopeMatrixAsStringList is not None
    assert extractor.updateIsotopeMatrixFromStringList is not None
    assert extractor.translateIsotopeMatrix is not None

@report
def testItraqQuantifier():
    """
    @tests:
     ItraqQuantifier.__init__
    """
    ff = pyopenms.ItraqQuantifier()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.ItraqQuantifier().run is not None

    # Note that using TMT_SIXPLEX will not work here
    p = ff.getDefaults()
    pyopenms.ItraqQuantifier(pyopenms.ITRAQ_TYPES.FOURPLEX, p)
    assert pyopenms.ItraqChannelExtractor().run is not None

    assert ff.getIsotopeMatrixAsStringList is not None
    assert ff.updateIsotopeMatrixFromStringList is not None
    assert ff.translateIsotopeMatrix is not None

@report
def testLinearResampler():
    """
    @tests:
     LinearResampler.__init__
    """
    ff = pyopenms.LinearResampler()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.LinearResampler().raster is not None
    assert pyopenms.LinearResampler().rasterExperiment is not None

@report
def testPeptideAndProteinQuant():
    """
    @tests:
     PeptideAndProteinQuant.__init__
    """
    ff = pyopenms.PeptideAndProteinQuant()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.PeptideAndProteinQuant().quantifyPeptides is not None
    assert pyopenms.PeptideAndProteinQuant().quantifyProteins is not None

@report
def testSeedListGenerator():
    """
    @tests:
     SeedListGenerator.__init__
    """
    ff = pyopenms.SeedListGenerator()
    p = ff.getDefaults()
    _testParam(p)

    # TODO 
    # assert pyopenms.SeedListGenerator().compute is not None

@report
def testTOFCalibration():
    """
    @tests:
     TOFCalibration.__init__
    """
    ff = pyopenms.TOFCalibration()
    p = ff.getDefaults()
    # _testParam(p)

    assert pyopenms.TOFCalibration().calibrate is not None
    assert pyopenms.TOFCalibration().pickAndCalibrate is not None

@report
def testConsensusID():
    """
    @tests:
     ConsensusID.__init__
    """
    ff = pyopenms.ConsensusID()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.ConsensusID().apply is not None

@report
def testFalseDiscoveryRate():
    """
    @tests:
     ConsensusID.__init__
    """
    ff = pyopenms.FalseDiscoveryRate()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.FalseDiscoveryRate().apply is not None

@report
def testIDFilter():
    """
    @tests:
     IDFilter.__init__
    """
    ff = pyopenms.IDFilter()

    # assert pyopenms.IDFilter().apply is not None

@report
def testProteinResolver():
    """
    @tests:
     ProteinResolver.__init__
    """
    ff = pyopenms.ProteinResolver()

    assert pyopenms.ProteinResolver().resolveConsensus is not None
    assert pyopenms.ProteinResolver().resolveID is not None
    assert pyopenms.ProteinResolver().setProteinData is not None
    assert pyopenms.ProteinResolver().getResults is not None

@report
def testSvmTheoreticalSpectrumGeneratorTrainer():
    """
    @tests:
     SvmTheoreticalSpectrumGeneratorTrainer.__init__
    """
    ff = pyopenms.SvmTheoreticalSpectrumGeneratorTrainer()

    assert pyopenms.SvmTheoreticalSpectrumGeneratorTrainer().trainModel is not None
    assert pyopenms.SvmTheoreticalSpectrumGeneratorTrainer().normalizeIntensity is not None

@report
def testPosteriorErrorProbabilityModel():
    """
    @tests:
     PosteriorErrorProbabilityModel.__init__
    """
    model = pyopenms.PosteriorErrorProbabilityModel()
    p = model.getDefaults()
    _testParam(p)

    assert pyopenms.PosteriorErrorProbabilityModel().fit is not None
    assert pyopenms.PosteriorErrorProbabilityModel().computeProbability is not None

    scores = [float(i) for i in range(10)]
    model.fit(scores)
    model.fit(scores, scores)

    model.fillDensities(scores, scores, scores)

    assert model.computeMaxLikelihood is not None
    assert model.one_minus_sum_post is not None
    assert model.sum_post is not None
    assert model.sum_pos_x0 is not None
    assert model.sum_neg_x0 is not None
    assert model.sum_pos_sigma is not None
    assert model.sum_neg_sigma is not None

    GaussFitResult = model.getCorrectlyAssignedFitResult()
    GaussFitResult = model.getIncorrectlyAssignedFitResult()
    model.getNegativePrior()
    model.getGauss(5.0, GaussFitResult)
    model.getGumbel(5.0, GaussFitResult)
    model.computeProbability(5.0) 

    # model.InitPlots

    target = [float(i) for i in range(10)]
    model.getGumbelGnuplotFormula(GaussFitResult) 
    model.getGaussGnuplotFormula(GaussFitResult) 
    model.getBothGnuplotFormula(GaussFitResult, GaussFitResult) 
    model.plotTargetDecoyEstimation(target, target)
    model.getSmallestScore()

@report
def testSeedListGenerator():
    """
    @tests:
     SeedListGenerator.__init__
    """
    ff = pyopenms.SeedListGenerator()

    # TODO 
    # assert pyopenms.SeedListGenerator().generateSeedList is not None

@report
def testConsensusMapNormalizerAlgorithmMedian():
    """
    @tests:
     ConsensusMapNormalizerAlgorithmMedian.__init__
    """
    ff = pyopenms.ConsensusMapNormalizerAlgorithmMedian()

    assert pyopenms.ConsensusMapNormalizerAlgorithmMedian().computeNormalizationFactors is not None
    assert pyopenms.ConsensusMapNormalizerAlgorithmMedian().normalizeMaps is not None

@report
def testConsensusMapNormalizerAlgorithmQuantile():
    """
    @tests:
     ConsensusMapNormalizerAlgorithmQuantile.__init__
    """
    ff = pyopenms.ConsensusMapNormalizerAlgorithmQuantile()

    assert pyopenms.ConsensusMapNormalizerAlgorithmQuantile().normalizeMaps is not None

@report
def testConsensusMapNormalizerAlgorithmThreshold():
    """
    @tests:
     ConsensusMapNormalizerAlgorithmThreshold.__init__
    """
    ff = pyopenms.ConsensusMapNormalizerAlgorithmThreshold()

    assert pyopenms.ConsensusMapNormalizerAlgorithmThreshold().computeCorrelation is not None
    assert pyopenms.ConsensusMapNormalizerAlgorithmThreshold().normalizeMaps is not None


@report
def testFeatureFinderAlgorithmPicked():
    """
    @tests:
     FeatureFinderAlgorithmPicked.__init__
    """
    ff = pyopenms.FeatureFinderAlgorithmPicked()

    assert pyopenms.FeatureFinderAlgorithmPicked().setData is not None
    assert pyopenms.FeatureFinderAlgorithmPicked().run is not None

@report
def testFeatureFinderAlgorithmSH():
    """
    @tests:
     FeatureFinderAlgorithmSH.__init__
    """
    ff = pyopenms.FeatureFinderAlgorithmSH()

    assert pyopenms.FeatureFinderAlgorithmSH().setData is not None
    assert pyopenms.FeatureFinderAlgorithmSH().run is not None

@report
def testFeatureFinderAlgorithmIsotopeWavelet():
    """
    @tests:
     FeatureFinderAlgorithmIsotopeWavelet.__init__
    """
    ff = pyopenms.FeatureFinderAlgorithmIsotopeWavelet()

    assert pyopenms.FeatureFinderAlgorithmIsotopeWavelet().setData is not None
    assert pyopenms.FeatureFinderAlgorithmIsotopeWavelet().run is not None


@report
def testAScore():
    """
    @tests:
     AScore.__init__
    """
    ff = pyopenms.AScore()

    hit = pyopenms.PeptideHit()
    richspectrum = pyopenms.RichMSSpectrum()

    ff.compute(hit, richspectrum, 5.0, 1)
    ff.computeCumulativeScore(1,1,0.5)

@report
def testIDRipper():
    """
    @tests:
     IDRipper.__init__
     IDRipper.rip
    """
    ff = pyopenms.IDRipper()

    assert pyopenms.IDRipper().rip is not None

@report
def testFASTAFile():
    """
    @tests:
     FASTAFile.__init__
     FASTAFile.load
     FASTAFile.store
    """
    ff = pyopenms.FASTAFile()

    assert pyopenms.FASTAFile().load is not None
    assert pyopenms.FASTAFile().store is not None


@report
def testFASTAEntry():
    """
    @tests:
     FASTAEntry.__init__
    """
    ff = pyopenms.FASTAEntry()

@report
def testInternalCalibration():
    """
    @tests:
     InternalCalibration.__init__
     InternalCalibration.calibrateMapGlobally
     InternalCalibration.calibrateMapSpectrumwise
    """
    ff = pyopenms.InternalCalibration()

    assert pyopenms.InternalCalibration().calibrateMapSpectrumwise is not None
    assert pyopenms.InternalCalibration().calibrateMapGlobally is not None

@report
def testTransitionTSVReader():
    """
    @tests:
     TransitionTSVReader.__init__
     TransitionTSVReader.calibrateMapGlobally
     TransitionTSVReader.calibrateMapSpectrumwise
    """
    ff = pyopenms.TransitionTSVReader()

    assert pyopenms.TransitionTSVReader().convertTargetedExperimentToTSV is not None
    assert pyopenms.TransitionTSVReader().convertTSVToTargetedExperiment is not None
    assert pyopenms.TransitionTSVReader().validateTargetedExperiment is not None

@report
def testEnzymaticDigestion():
    """
    @tests:
     EnzymaticDigestion.__init__
     EnzymaticDigestion.getMissedCleavages()
     EnzymaticDigestion.setMissedCleavages()
     EnzymaticDigestion.getEnzyme()
     EnzymaticDigestion.setEnzyme()
     EnzymaticDigestion.getEnzymeByName()
     EnzymaticDigestion.digest()
     EnzymaticDigestion.peptideCount()
     EnzymaticDigestion.isLogModelEnabled()
     EnzymaticDigestion.setLogModelEnabled()
     EnzymaticDigestion.getLogThreshold()
     EnzymaticDigestion.setLogThreshold()
    """
    ff = pyopenms.EnzymaticDigestion()
    enz = pyopenms.EnzymaticDigestion().Enzyme()

    assert pyopenms.EnzymaticDigestion().getMissedCleavages is not None
    assert pyopenms.EnzymaticDigestion().setMissedCleavages is not None
    assert pyopenms.EnzymaticDigestion().getEnzyme is not None
    assert pyopenms.EnzymaticDigestion().setEnzyme is not None
    assert pyopenms.EnzymaticDigestion().getEnzymeByName is not None

    assert pyopenms.EnzymaticDigestion().digest is not None
    assert pyopenms.EnzymaticDigestion().peptideCount is not None

    assert pyopenms.EnzymaticDigestion().isLogModelEnabled is not None
    assert pyopenms.EnzymaticDigestion().setLogModelEnabled is not None
    assert pyopenms.EnzymaticDigestion().getLogThreshold  is not None
    assert pyopenms.EnzymaticDigestion().setLogThreshold is not None

    ff.setLogThreshold(5) 
    assert ff.getLogThreshold() == 5 

    ff.setMissedCleavages(5) 
    assert ff.getMissedCleavages() == 5 

    ff.setEnzyme(enz.TRYPSIN) 
    assert ff.getEnzyme() == enz.TRYPSIN

    ff.setLogModelEnabled(True) 
    assert ff.isLogModelEnabled() == True

@report
def testIDDecoyProbability():
    """
    @tests:
      IDDecoyProbability.__init__
    """
    ff = pyopenms.IDDecoyProbability()

    assert pyopenms.IDDecoyProbability().apply is not None

@report
def testFeatureGrouping():
    """
    @tests:
     FeatureGroupingAlgorithm.getDefaults
     FeatureGroupingAlgorithm.getName
     FeatureGroupingAlgorithm.getParameters
     FeatureGroupingAlgorithm.setName
     FeatureGroupingAlgorithm.setParameters
     FeatureGroupingAlgorithm.transferSubelements
     FeatureGroupingAlgorithmQT.__init__
     FeatureGroupingAlgorithmQT.getDefaults
     FeatureGroupingAlgorithmQT.getName
     FeatureGroupingAlgorithmQT.getParameters
     FeatureGroupingAlgorithmQT.group
     FeatureGroupingAlgorithmQT.setName
     FeatureGroupingAlgorithmQT.setParameters
     FeatureGroupingAlgorithmQT.transferSubelements
    """

    assert pyopenms.FeatureGroupingAlgorithm.getDefaults is not None
    assert pyopenms.FeatureGroupingAlgorithm.getName is not None
    assert pyopenms.FeatureGroupingAlgorithm.getParameters is not None
    assert pyopenms.FeatureGroupingAlgorithm.setName is not None
    assert pyopenms.FeatureGroupingAlgorithm.setParameters is not None
    assert pyopenms.FeatureGroupingAlgorithm.transferSubelements is not None

    qt = pyopenms.FeatureGroupingAlgorithmQT()
    qt.getDefaults()
    qt.getParameters()
    qt.getName()
    assert qt.group is not None
    assert qt.setName is not None
    assert qt.setParameters is not None
    assert qt.transferSubelements is not None

@report
def testFeatureMap():
    """
    @tests:
     FeatureMap.__init__
     FeatureMap.__add__
     FeatureMap.__iadd__
     FeatureMap.__radd__
     FeatureMap.__getitem__
     FeatureMap.__iter__
     FeatureMap.clear
     FeatureMap.clearUniqueId
     FeatureMap.ensureUniqueId
     FeatureMap.getDataProcessing
     FeatureMap.getProteinIdentifications
     FeatureMap.getUnassignedPeptideIdentifications
     FeatureMap.getUniqueId
     FeatureMap.setUniqueId
     FeatureMap.hasInvalidUniqueId
     FeatureMap.hasValidUniqueId
     FeatureMap.push_back
     FeatureMap.setDataProcessing
     FeatureMap.setProteinIdentifications
     FeatureMap.setUnassignedPeptideIdentifications
     FeatureMap.setUniqueIds
     FeatureMap.size
     FeatureMap.sortByIntensity
     FeatureMap.sortByMZ
     FeatureMap.sortByOverallQuality
     FeatureMap.sortByPosition
     FeatureMap.sortByRT
     FeatureMap.swap
     FeatureMap.updateRanges
    """
    fm = pyopenms.FeatureMap()
    _testUniqueIdInterface(fm)
    fm.clear()
    fm.clearUniqueId()

    fm.getIdentifier()
    fm.getLoadedFileType()
    fm.getLoadedFilePath()

    f = pyopenms.Feature()
    fm.push_back(f)

    assert len(list(fm)) == 1


    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByIntensity()
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByIntensity(False)
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByPosition()
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByRT()
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByMZ()
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByOverallQuality()
    assert fm.size() == 1
    assert fm[0] == f

    fm2 = pyopenms.FeatureMap()

    fm.swap(fm2)
    assert fm2.size() == 1
    assert fm2[0] == f

    assert fm.size() == 0

    fm2.updateRanges()

    assert isinstance(fm2.getMin()[0], float)
    assert isinstance(fm2.getMin()[1], float)
    assert isinstance(fm2.getMax()[0], float)
    assert isinstance(fm2.getMax()[1], float)
    assert isinstance(fm2.getMinInt(), float)
    assert isinstance(fm2.getMaxInt(), float)

    assert fm2.getProteinIdentifications() == []
    fm2.setProteinIdentifications([])

    assert fm2.getUnassignedPeptideIdentifications() == []
    fm2.setUnassignedPeptideIdentifications([])

    fm2.clear()
    assert fm2.size() == 0

    dp = pyopenms.DataProcessing()
    fm2.setDataProcessing([dp])
    assert fm2.getDataProcessing() == [dp]
    testDataProcessing(dp)

    fm2.setUniqueIds()

    fm += fm
    assert fm + fm != fm


@report
def testFeatureXMLFile():
    """
    @tests:
     FeatureXMLFile.__init__
     FeatureXMLFile.load
     FeatureXMLFile.store
     FeatureXMLFile.getOptions
     FeatureXMLFile.setOptions
     FeatureXMLFile.loadSize

     FileHandler.__init__
     FileHandler.loadFeature
    """

    fm = pyopenms.FeatureMap()
    fm.setUniqueIds()
    fh = pyopenms.FeatureXMLFile()
    fh.store("test.featureXML", fm)
    fh.load("test.featureXML", fm)

    fh = pyopenms.FileHandler()
    fh.loadFeatures("test.featureXML", fm)

@report
def testFileDescription():
    """
    @tests:
     FileDescription.__init__
     FileDescription.filename
     FileDescription.label
     FileDescription.size
     FileDescription.unique_id
    """
    fd = pyopenms.FileDescription()
    assert isinstance(fd.filename, str)
    assert isinstance(fd.label, str)
    assert isinstance(fd.size, int)
    assert isinstance(fd.unique_id, (long, int, str))

@report
def testFileHandler():
    """
    @tests:
     FileHandler.__init__
     FileHandler.getType
     FileHandler.loadExperiment
     FileHandler.storeExperiment
    """
    mse = pyopenms.MSExperiment()

    fh = pyopenms.FileHandler()
    fh.storeExperiment("test1.mzML", mse)
    fh.loadExperiment("test1.mzML", mse)
    fh.storeExperiment("test1.mzXML", mse)
    fh.loadExperiment("test1.mzXML", mse)
    fh.storeExperiment("test1.mzData", mse)
    fh.loadExperiment("test1.mzData", mse)

    

@report
def testIDMapper():
    """
    @tests:
     IDMapper.__init__
     IDMapper.annotate
     IDMapper.getDefaults
     IDMapper.getName
     IDMapper.getParameters
     IDMapper.setName
     IDMapper.setParameters
    """
    idm = pyopenms.IDMapper()
    assert idm.annotate is not None
    idm.getDefaults()
    idm.setName("x")
    assert idm.getName() == "x"
    idm.setParameters(idm.getParameters())

@report
def testIdXMLFile():
    """
    @tests:
     IdXMLFile.__init__
     IdXMLFile.load
     IdXMLFile.store
    """
    assert pyopenms.IdXMLFile().load is not None
    assert pyopenms.IdXMLFile().store is not None

@report
def testPepXMLFile():
    """
    @tests:
     PepXMLFile.__init__
     PepXMLFile.load
     PepXMLFile.store
    """
    f = pyopenms.PepXMLFile()

    assert pyopenms.PepXMLFile().load is not None
    assert pyopenms.PepXMLFile().store is not None

@report
def testProtXMLFile():
    """
    @tests:
     ProtXMLFile.__init__
     ProtXMLFile.load
     ProtXMLFile.store
    """
    f = pyopenms.ProtXMLFile()

    assert pyopenms.ProtXMLFile().load is not None
    assert pyopenms.ProtXMLFile().store is not None

@report
def testDTA2DFile():
    """
    @tests:
     DTA2DFile.__init__
     DTA2DFile.load
     DTA2DFile.store
    """
    f = pyopenms.DTA2DFile()

    assert pyopenms.DTA2DFile().load is not None
    assert pyopenms.DTA2DFile().store is not None

@report
def testDTAFile():
    """
    @tests:
     DTAFile.__init__
     DTAFile.load
     DTAFile.store
    """
    f = pyopenms.DTAFile()

    assert pyopenms.DTAFile().load is not None
    assert pyopenms.DTAFile().store is not None

@report
def testEDTAFile():
    """
    @tests:
     EDTAFile.__init__
     EDTAFile.load
     EDTAFile.store
    """
    f = pyopenms.EDTAFile()

    assert pyopenms.EDTAFile().load is not None
    assert pyopenms.EDTAFile().store is not None

@report
def testKroenikFile():
    """
    @tests:
     KroenikFile.__init__
     KroenikFile.load
     KroenikFile.store
    """
    f = pyopenms.KroenikFile()

    assert pyopenms.KroenikFile().load is not None
    assert pyopenms.KroenikFile().store is not None

@report
def testMSPFile():
    """
    @tests:
     MSPFile.__init__
    """
    f = pyopenms.MSPFile()

    # assert pyopenms.KroenikFile().load is not None
    # assert pyopenms.KroenikFile().store is not None

@report
def testMzIdentMLFile():
    """
    @tests:
     MzIdentMLFile.__init__
    """
    f = pyopenms.MzIdentMLFile()

    assert pyopenms.MzIdentMLFile().load is not None
    assert pyopenms.MzIdentMLFile().store is not None
    assert pyopenms.MzIdentMLFile().isSemanticallyValid is not None


@report
def testMzTabFile():
    """
    @tests:
     MzTabFile.__init__
    """
    f = pyopenms.MzTabFile()

    # assert pyopenms.MzTabFile().store is not None

@report
def testMzTab():
    """
    @tests:
     MzTab.__init__
    """
    # f = pyopenms.MzTab()

@report
def testInstrumentSettings():
    """
    @tests:
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

@report
def testContactPerson():
    """
    @tests:
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
    @tests:
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
    @tests:
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
    @tests:
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
    @tests:
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

@report
def testIonDetector():
    """
    @tests:
     IonDetector.__init__
     IonDetector.setAcquisitionMode
     IonDetector.getAcquisitionMode
     IonDetector.setResolution
     IonDetector.getResolution
     IonDetector.setADCSamplingFrequency
     IonDetector.getADCSamplingFrequency
     IonDetector.setOrder
     IonDetector.getOrder
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

@report
def testIonSource():
    """
    @tests:
     IonSource.__init__
     IonSource.setPolarity
     IonSource.getPolarity
     IonSource.setInletType
     IonSource.getInletType
     IonSource.setIonizationMethod
     IonSource.getIonizationMethod
     IonSource.setOrder
     IonSource.getOrder
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

@report
def testMassAnalyzer():
    """
    @tests:
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

@report
def testSample():
    """
    @tests:
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
     Sample.removeTreatment
     Sample.countTreatments
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

    has_exception = False
    try:
        ins.removeTreatment(0)
    except Exception:
        has_exception = True
    assert has_exception
    assert ins.countTreatments() == 0

@report
def testLogType():

    """
    @tests:
     LogType.CMD
     LogType.GUI
     LogType.NONE
     """
    assert isinstance(pyopenms.LogType.CMD, int)
    assert isinstance(pyopenms.LogType.GUI, int)
    assert isinstance(pyopenms.LogType.NONE, int)

@report
def testMSExperiment():
    """
    @tests:
     MSExperiment.__init__
     MSExperiment.getLoadedFilePath
     MSExperiment.getMaxMZ
     MSExperiment.getMaxRT
     MSExperiment.getMetaValue
     MSExperiment.getMinMZ
     MSExperiment.getMinRT
     MSExperiment.push_back
     MSExperiment.setLoadedFilePath
     MSExperiment.setMetaValue
     MSExperiment.size
     MSExperiment.sortSpectra
     MSExperiment.updateRanges
     MSExperiment.__eq__
     MSExperiment.__ge__
     MSExperiment.__getitem__
     MSExperiment.__gt__
     MSExperiment.__iter__
     MSExperiment.__le__
     MSExperiment.__lt__
     MSExperiment.__ne__
     MSExperiment.clearMetaInfo
     MSExperiment.getKeys
     MSExperiment.isMetaEmpty
     MSExperiment.metaValueExists
     MSExperiment.removeMetaValue
     MSExperiment.getSize
     MSExperiment.isSorted
    """
    mse = pyopenms.MSExperiment()
    _testMetaInfoInterface(mse)
    mse.updateRanges()
    mse.sortSpectra(True)
    assert isinstance(mse.getMaxRT(), float)
    assert isinstance(mse.getMinRT(), float)
    assert isinstance(mse.getMaxMZ(), float)
    assert isinstance(mse.getMinMZ(), float)
    assert isinstance(mse.getLoadedFilePath(), str)
    assert isinstance(mse.getMinInt(), float)
    assert isinstance(mse.getMaxInt(), float)

    assert isinstance(mse.getMin()[0], float)
    assert isinstance(mse.getMin()[1], float)
    assert isinstance(mse.getMax()[0], float)
    assert isinstance(mse.getMax()[1], float)
    mse.setLoadedFilePath("")
    assert mse.size() == 0

    mse.getIdentifier()
    mse.getLoadedFileType()
    mse.getLoadedFilePath()

    mse.addSpectrum(pyopenms.MSSpectrum())
    assert mse.size() == 1

    assert mse[0] is not None

    assert isinstance(list(mse), list)

    assert mse == mse
    assert not mse != mse

    assert mse.getSize() >= 0
    assert int(mse.isSorted()) in (0,1)

    import copy
    mse2 = copy.copy(mse)

    assert mse.getSize() == mse2.getSize()
    assert mse2 == mse


@report
def testMSQuantifications():
    """
    @tests:
     MSQuantifications.__eq__
     MSQuantifications.__ge__
     MSQuantifications.__gt__
     MSQuantifications.__init__
     MSQuantifications.__le__
     MSQuantifications.__lt__
     MSQuantifications.__ne__
     MSQuantifications.getConsensusMaps
     MSQuantifications.setConsensusMaps
     MSQuantifications.setDataProcessing
     MSQuantifications.getDataProcessing
     MSQuantifications.getAssays
     MSQuantifications.getFeatureMaps
     MSQuantifications.setAnalysisSummaryQuantType
     MSQuantifications.getAnalysisSummary
     MSQuantifications.addConsensusMap
     MSQuantifications.assignUIDs
    """
    msq = pyopenms.MSQuantifications()
    assert msq == msq
    assert not msq != msq
    msq.setConsensusMaps(msq.getConsensusMaps())

    summary = msq.getAnalysisSummary()
    msq.setDataProcessingList(msq.getDataProcessingList())
    msq.getAssays()
    msq.getFeatureMaps()
    msq.setAnalysisSummaryQuantType(pyopenms.MSQuantifications.QUANT_TYPES.LABELFREE)

    msq.addConsensusMap(pyopenms.ConsensusMap())
    msq.assignUIDs()

@report
def testMSSpectrum():
    """
    @tests:
     MSSpectrum.__init__
     MSSpectrum.clear
     MSSpectrum.clearMetaInfo
     MSSpectrum.findNearest
     MSSpectrum.getAcquisitionInfo
     MSSpectrum.getComment
     MSSpectrum.getDataProcessing
     MSSpectrum.getInstrumentSettings
     MSSpectrum.getKeys
     MSSpectrum.getMSLevel
     MSSpectrum.getMetaValue
     MSSpectrum.getName
     MSSpectrum.getNativeID
     MSSpectrum.getPeptideIdentifications
     MSSpectrum.getPrecursors
     MSSpectrum.getProducts
     MSSpectrum.getRT
     MSSpectrum.getSourceFile
     MSSpectrum.getType
     MSSpectrum.get_peaks
     MSSpectrum.intensityInRange
     MSSpectrum.isMetaEmpty
     MSSpectrum.isSorted
     MSSpectrum.metaValueExists
     MSSpectrum.push_back
     MSSpectrum.removeMetaValue
     MSSpectrum.setAcquisitionInfo
     MSSpectrum.setComment
     MSSpectrum.setDataProcessing
     MSSpectrum.setInstrumentSettings
     MSSpectrum.setMSLevel
     MSSpectrum.setMetaValue
     MSSpectrum.setName
     MSSpectrum.setNativeID
     MSSpectrum.setPeptideIdentifications
     MSSpectrum.setPrecursors
     MSSpectrum.setProducts
     MSSpectrum.setRT
     MSSpectrum.setSourceFile
     MSSpectrum.setType
     MSSpectrum.set_peaks
     MSSpectrum.size
     MSSpectrum.unify
     MSSpectrum.updateRanges
     MSSpectrum.__eq__
     MSSpectrum.__ge__
     MSSpectrum.__getitem__
     MSSpectrum.__gt__
     MSSpectrum.__le__
     MSSpectrum.__lt__
     MSSpectrum.__ne__
     """
    spec = pyopenms.MSSpectrum()
    _testMetaInfoInterface(spec)

    testSpectrumSetting(spec)

    spec.setRT(3.0)
    assert spec.getRT() == 3.0
    spec.setMSLevel(2)
    assert spec.getMSLevel() == 2
    spec.setName("spec")
    assert spec.getName() == "spec"

    p = pyopenms.Peak1D()
    p.setMZ(1000.0)
    p.setIntensity(200.0)

    spec.push_back(p)
    assert spec.size() == 1
    assert spec[0] == p

    spec.updateRanges()
    assert isinstance(spec.findNearest(0.0), int)

    assert isinstance(spec.getMin()[0], float)
    assert isinstance(spec.getMax()[0], float)
    assert isinstance(spec.getMinInt(), float)
    assert isinstance(spec.getMaxInt(), float)

    assert spec == spec
    assert not spec != spec

    assert spec.get_peaks().shape == (1,2), spec.get_peaks().shape

    assert int(spec.isSorted()) in  (0,1)

@report
def testMRMFeature():
    """
    @tests:
      MRMFeature.__init__
      MRMFeature.addScore
      MRMFeature.getScore
     """
    mrmfeature = pyopenms.MRMFeature()

    mrmfeature.addScore("testscore", 6)
    assert mrmfeature.getScore("testscore") == 6.0
    mrmfeature.addScore("testscore", 7)
    assert mrmfeature.getScore("testscore") == 7.0

@report
def testConfidenceScoring():
    """
    @tests:
      ConfidenceScoring.__init__
     """
    scoring = pyopenms.ConfidenceScoring()

@report
def testMRMDecoy():
    """
    @tests:
      MRMDecoy.__init__
     """
    mrmdecoy = pyopenms.MRMDecoy()
    assert mrmdecoy is not None

    assert pyopenms.MRMDecoy().restrictTransitions is not None
    assert pyopenms.MRMDecoy().generateDecoys is not None

@report
def testMRMTransitionGroup():
    """
    @tests:
     """
    mrmgroup = pyopenms.MRMTransitionGroup()
    assert mrmgroup is not None

    mrmgroup.setTransitionGroupID("this_id")
    assert mrmgroup.getTransitionGroupID() == "this_id"

    assert len(mrmgroup.getTransitions()) == 0
    mrmgroup.addTransition(pyopenms.ReactionMonitoringTransition(), "tr1")
    assert len(mrmgroup.getTransitions()) == 1

@report
def testReactionMonitoringTransition():
    """
    @tests:
     """
    tr = pyopenms.ReactionMonitoringTransition()

@report
def testMapAlignment():

    """
    @tests:
     MapAlignmentAlgorithmPoseClustering.__init__
     MapAlignmentAlgorithmPoseClustering.getDefaults
     MapAlignmentAlgorithmPoseClustering.getName
     MapAlignmentAlgorithmPoseClustering.getParameters
     MapAlignmentAlgorithmPoseClustering.setName
     MapAlignmentAlgorithmPoseClustering.setParameters
     MapAlignmentAlgorithmPoseClustering.setReference

     MapAlignmentAlgorithmPoseClustering.align
     MapAlignmentAlgorithmPoseClustering.endProgress
     MapAlignmentAlgorithmPoseClustering.getLogType
     MapAlignmentAlgorithmPoseClustering.setLogType
     MapAlignmentAlgorithmPoseClustering.setProgress
     MapAlignmentAlgorithmPoseClustering.startProgress

     MapAlignmentTransformer.transformFeatureMaps
     MapAlignmentTransformer.transformPeakMaps
     MapAlignmentTransformer.transformSingleFeatureMap
     MapAlignmentTransformer.transformSinglePeakMap
     """
    ma = pyopenms.MapAlignmentAlgorithmPoseClustering()
    assert isinstance(ma.getDefaults(), pyopenms.Param)
    assert isinstance(ma.getParameters(), pyopenms.Param)
    assert isinstance(ma.getName(), str)

    ma.setName(ma.getName())

    ma.getDefaults()
    ma.getParameters()

    ma.setParameters(ma.getDefaults())

    ma.setReference
    ma.align


    pyopenms.MapAlignmentTransformer.transformPeakMaps
    pyopenms.MapAlignmentTransformer.transformFeatureMaps
    pyopenms.MapAlignmentTransformer.transformSinglePeakMap
    pyopenms.MapAlignmentTransformer.transformSingleFeatureMap

@report
def testMapAlignmentIdentification():

    """
    @tests:
     MapAlignmentAlgorithmIdentification.__init__
     """
    ma = pyopenms.MapAlignmentAlgorithmIdentification()

    assert pyopenms.MapAlignmentAlgorithmIdentification().alignPeakMaps is not None
    assert pyopenms.MapAlignmentAlgorithmIdentification().alignFeatureMaps is not None
    assert pyopenms.MapAlignmentAlgorithmIdentification().alignConsensusMaps is not None
    assert pyopenms.MapAlignmentAlgorithmIdentification().alignPeptideIdentifications is not None
    assert pyopenms.MapAlignmentAlgorithmIdentification().setReference is not None
    assert pyopenms.MapAlignmentAlgorithmIdentification().fitModel is not None

@report
def testMapAlignmentTransformer():

    """
    @tests:
     MapAlignmentTransformer.__init__
     """
    ma = pyopenms.MapAlignmentTransformer()

    assert pyopenms.MapAlignmentTransformer().transformPeakMaps is not None
    assert pyopenms.MapAlignmentTransformer().transformFeatureMaps is not None
    assert pyopenms.MapAlignmentTransformer().transformConsensusMaps is not None
    # assert pyopenms.MapAlignmentTransformer().transformPeptideIdentifications is not None
    assert pyopenms.MapAlignmentTransformer().transformSinglePeakMap is not None
    assert pyopenms.MapAlignmentTransformer().transformSingleFeatureMap is not None
    assert pyopenms.MapAlignmentTransformer().transformSingleConsensusMap is not None
    assert pyopenms.MapAlignmentTransformer().transformSinglePeptideIdentification is not None


@report
def testMxxxFile():
    """
    @tests:
     MzDataFile.__init__
     MzDataFile.endProgress
     MzDataFile.getLogType
     MzDataFile.load
     MzDataFile.setLogType
     MzDataFile.setProgress
     MzDataFile.startProgress
     MzDataFile.store
     MzDataFile.getOptions
     MzDataFile.setOptions

     MzMLFile.__init__
     MzMLFile.endProgress
     MzMLFile.getLogType
     MzMLFile.load
     MzMLFile.setLogType
     MzMLFile.setProgress
     MzMLFile.startProgress
     MzMLFile.store
     MzMLFile.getOptions
     MzMLFile.setOptions

     MzXMLFile.getOptions
     MzXMLFile.setOptions
     MzXMLFile.__init__
     MzXMLFile.endProgress
     MzXMLFile.getLogType
     MzXMLFile.load
     MzXMLFile.setLogType
     MzXMLFile.setProgress
     MzXMLFile.startProgress
     MzXMLFile.store

     MzQuantMLFile.__init__
     MzQuantMLFile.isSemanticallyValid
     MzQuantMLFile.load
     MzQuantMLFile.store
    """
    mse = pyopenms.MSExperiment()

    fh = pyopenms.MzDataFile()
    _testProgressLogger(fh)
    fh.store("test.mzData", mse)
    fh.load("test.mzData", mse)

    fh.setOptions(fh.getOptions())

    fh = pyopenms.MzMLFile()
    _testProgressLogger(fh)
    fh.store("test.mzML", mse)
    fh.load("test.mzML", mse)
    fh.setOptions(fh.getOptions())

    fh = pyopenms.MzXMLFile()
    _testProgressLogger(fh)
    fh.store("test.mzXML", mse)
    fh.load("test.mzXML", mse)
    fh.setOptions(fh.getOptions())

    fh = pyopenms.MzQuantMLFile()
    fh.isSemanticallyValid
    fh.load
    fh.store



@report
def testParamXMLFile():

    """
    @tests:
     ParamXMLFile.__init__
     ParamXMLFile.load
     ParamXMLFile.store
    """

    fh = pyopenms.ParamXMLFile()
    p = pyopenms.Param()
    fh.store("test.ini", p)
    fh.load("test.ini", p)



@report
def testPeak():

    """
    @tests:
     Peak1D.__init__
     Peak1D.getIntensity
     Peak1D.getMZ
     Peak1D.setIntensity
     Peak1D.setMZ
     Peak1D.__eq__
     Peak1D.__ge__
     Peak1D.__gt__
     Peak1D.__le__
     Peak1D.__lt__
     Peak1D.__ne__
     Peak2D.__init__
     Peak2D.getIntensity
     Peak2D.getMZ
     Peak2D.getRT
     Peak2D.setIntensity
     Peak2D.setMZ
     Peak2D.setRT
     Peak2D.__eq__
     Peak2D.__ge__
     Peak2D.__gt__
     Peak2D.__le__
     Peak2D.__lt__
     Peak2D.__ne__
    """
    p1 = pyopenms.Peak1D()
    p1.setIntensity(12.0)
    assert p1.getIntensity() == 12.0
    p1.setMZ(13.0)
    assert p1.getMZ() == 13.0

    assert p1 == p1
    assert not p1 != p1

    p2 = pyopenms.Peak2D()
    assert p2 == p2
    assert not p2 != p2
    p2.setIntensity(22.0)
    assert p2.getIntensity() == 22.0
    p2.setMZ(23.0)
    assert p2.getMZ() == 23.0
    p2.setRT(45.0)
    assert p2.getRT() == 45.0



@report
def testPeakFileOptions():
    """
    @tests:
     PeakFileOptions.__init__
     PeakFileOptions.addMSLevel
     PeakFileOptions.clearMSLevels
     PeakFileOptions.containsMSLevel
     PeakFileOptions.getCompression
     PeakFileOptions.getMSLevels
     PeakFileOptions.getMetadataOnly
     PeakFileOptions.getWriteSupplementalData
     PeakFileOptions.hasMSLevels
     PeakFileOptions.setCompression
     PeakFileOptions.setMSLevels
     PeakFileOptions.setMetadataOnly
     PeakFileOptions.setWriteSupplementalData
    """

    pfo = pyopenms.PeakFileOptions()
    pfo.addMSLevel
    pfo.clearMSLevels()
    pfo.containsMSLevel(1)
    pfo.getCompression()
    pfo.getMSLevels()
    pfo.getMetadataOnly()
    pfo.getWriteSupplementalData()
    pfo.hasMSLevels()
    pfo.setCompression
    pfo.setMSLevels
    pfo.setMetadataOnly
    pfo.setWriteSupplementalData


@report
def testPeakPickerHiRes():
    """
    @tests:
     PeakPickerHiRes.__init__
     PeakPickerHiRes.endProgress
     PeakPickerHiRes.getDefaults
     PeakPickerHiRes.getLogType
     PeakPickerHiRes.getName
     PeakPickerHiRes.getParameters
     PeakPickerHiRes.pick
     PeakPickerHiRes.pickExperiment
     PeakPickerHiRes.setLogType
     PeakPickerHiRes.setName
     PeakPickerHiRes.setParameters
     PeakPickerHiRes.setProgress
     PeakPickerHiRes.startProgress
    """

@report
def testPeakTypeEstimator():
    """
    @tests:
     PeakTypeEstimator.__init__
     PeakTypeEstimator.estimateType
    """

    pyopenms.PeakTypeEstimator().estimateType(pyopenms.MSSpectrum())

@report
def testPeptideHit():
    """
    @tests:
     PeptideHit.__init__
     PeptideHit.addProteinAccession
     PeptideHit.clearMetaInfo
     PeptideHit.getAAAfter
     PeptideHit.getAABefore
     PeptideHit.getKeys
     PeptideHit.getMetaValue
     PeptideHit.getProteinAccessions
     PeptideHit.getRank
     PeptideHit.getScore
     PeptideHit.getSequence
     PeptideHit.isMetaEmpty
     PeptideHit.metaValueExists
     PeptideHit.removeMetaValue
     PeptideHit.setAAAfter
     PeptideHit.setAABefore
     PeptideHit.setCharge
     PeptideHit.setMetaValue
     PeptideHit.setProteinAccessions
     PeptideHit.setRank
     PeptideHit.setScore
     PeptideHit.setSequence
     PeptideHit.__eq__
     PeptideHit.__ge__
     PeptideHit.__gt__
     PeptideHit.__le__
     PeptideHit.__lt__
     PeptideHit.__ne__
    """
    ph = pyopenms.PeptideHit()
    assert ph == ph
    assert not ph != ph

    ph = pyopenms.PeptideHit(1.0, 1, 0, pyopenms.AASequence("A"))
    _testMetaInfoInterface(ph)
    ph.addProteinAccession("A")
    assert ph.getProteinAccessions() == ["A"]

    assert ph.getScore() == 1.0
    assert ph.getRank() == 1
    assert ph.getSequence().toString() == "A"

    ph.setScore(2.0)
    assert ph.getScore() == 2.0
    ph.setRank(30)
    assert ph.getRank() == 30
    ph.setSequence(pyopenms.AASequence("AAA"))
    assert ph.getSequence().toString() == "AAA"

    ph.setAABefore('B')
    assert ph.getAABefore() == "B"
    ph.setAAAfter('C')
    assert ph.getAAAfter() == 'C'

    assert ph == ph
    assert not ph != ph


@report
def testPeptideIdentification():
    """
    @tests:
     PeptideIdentification.__init__
     PeptideIdentification.assignRanks
     PeptideIdentification.clearMetaInfo
     PeptideIdentification.empty
     PeptideIdentification.getHits
     PeptideIdentification.getIdentifier
     PeptideIdentification.getKeys
     PeptideIdentification.getMetaValue
     PeptideIdentification.getNonReferencingHits
     PeptideIdentification.getReferencingHits
     PeptideIdentification.getScoreType
     PeptideIdentification.getSignificanceThreshold
     PeptideIdentification.insertHit
     PeptideIdentification.isHigherScoreBetter
     PeptideIdentification.isMetaEmpty
     PeptideIdentification.metaValueExists
     PeptideIdentification.removeMetaValue
     PeptideIdentification.setHigherScoreBetter
     PeptideIdentification.setHits
     PeptideIdentification.setIdentifier
     PeptideIdentification.setMetaValue
     PeptideIdentification.setScoreType
     PeptideIdentification.sort
     PeptideIdentification.__eq__
     PeptideIdentification.__ge__
     PeptideIdentification.__gt__
     PeptideIdentification.__le__
     PeptideIdentification.__lt__
     PeptideIdentification.__ne__
     PeptideIdentification.setSignificanceThreshold
     """
    pi = pyopenms.PeptideIdentification()
    _testMetaInfoInterface(pi)
    assert pi == pi
    assert not pi != pi

    ph = pyopenms.PeptideHit(1.0, 1, 0, pyopenms.AASequence("A"))
    pi.insertHit(ph)
    phx, = pi.getHits()
    assert phx == ph

    pi.setHits([ph])
    phx, = pi.getHits()
    assert phx == ph

    assert isinstance(pi.getSignificanceThreshold(), float)
    assert isinstance(pi.getScoreType(), str)
    pi.setScoreType("A")
    assert isinstance(pi.isHigherScoreBetter(), int)
    assert isinstance(pi.getIdentifier(), str)
    pi.setIdentifier("id")
    pi.assignRanks()
    pi.sort()
    assert not pi.empty()

    rv = []
    pi.getReferencingHits("A", rv)
    assert rv == []
    pi.getNonReferencingHits("A", rv)
    hit, = rv
    assert hit.getSequence().toString()== "A"
    assert hit.getScore() == 1.0
    assert hit.getRank() == 1

    rv = []
    pi.getReferencingHits(["A"], rv)
    assert rv == []
    pi.getNonReferencingHits(["A"], rv)
    hit, = rv
    assert hit.getSequence().toString()== "A"
    assert hit.getScore() == 1.0
    assert hit.getRank() == 1

    ph = pyopenms.ProteinHit()
    pi.getReferencingHits([ph], rv)
    hit, = rv
    assert hit.getSequence().toString()== "A"
    assert hit.getScore() == 1.0
    assert hit.getRank() == 1
    rv = []
    pi.getNonReferencingHits([ph], rv)
    hit, = rv
    assert hit.getSequence().toString()== "A"
    assert hit.getScore() == 1.0
    assert hit.getRank() == 1

    pi.setSignificanceThreshold(6.0)


@report
def testPolarity():
    """
    @tests:
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
    @tests:
     Precursor.__init__
     Precursor.getIntensity
     Precursor.getMZ
     Precursor.setIntensity
     Precursor.setMZ
     Precursor.setActivationMethods
     Precursor.getActivationMethods
     Precursor.setActivationEnergy
     Precursor.getActivationEnergy
     Precursor.setIsolationWindowUpperOffset
     Precursor.getIsolationWindowUpperOffset
     Precursor.setIsolationWindowLowerOffset
     Precursor.getIsolationWindowLowerOffset
     Precursor.setCharge
     Precursor.getCharge
     Precursor.setPossibleChargeStates
     Precursor.getPossibleChargeStates
     Precursor.getUnchargedMass
    """
    pc = pyopenms.Precursor()
    pc.setMZ(123.0)
    pc.setIntensity(12.0)
    assert pc.getMZ() == 123.0
    assert pc.getIntensity() == 12.0

    pc.setActivationMethods(pc.getActivationMethods())
    pc.setActivationEnergy(6.0)
    pc.getActivationEnergy()

    pc.setIsolationWindowUpperOffset(500.0)
    pc.getIsolationWindowUpperOffset()
    pc.setIsolationWindowLowerOffset(600.0)
    pc.getIsolationWindowLowerOffset()

    pc.setCharge(2)
    pc.getCharge()

    pc.setPossibleChargeStates(pc.getPossibleChargeStates())

    pc.getUnchargedMass()

@report
def testProcessingAction():
    """
    @tests:
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
     ProcessingAction.IDENTIFICATION_MAPPING
     ProcessingAction.NORMALIZATION
     ProcessingAction.PEAK_PICKING
     ProcessingAction.PRECURSOR_RECALCULATION
     ProcessingAction.QUANTITATION
     ProcessingAction.SIZE_OF_PROCESSINGACTION
     ProcessingAction.SMOOTHING
    """
    assert isinstance(pyopenms.ProcessingAction.ALIGNMENT, int)
    assert isinstance(pyopenms.ProcessingAction.BASELINE_REDUCTION, int)
    assert isinstance(pyopenms.ProcessingAction.CALIBRATION, int)
    assert isinstance(pyopenms.ProcessingAction.CHARGE_CALCULATION, int)
    assert isinstance(pyopenms.ProcessingAction.CHARGE_DECONVOLUTION, int)
    assert isinstance(pyopenms.ProcessingAction.CONVERSION_DTA, int)
    assert isinstance(pyopenms.ProcessingAction.CONVERSION_MZDATA, int)
    assert isinstance(pyopenms.ProcessingAction.CONVERSION_MZML, int)
    assert isinstance(pyopenms.ProcessingAction.CONVERSION_MZXML, int)
    assert isinstance(pyopenms.ProcessingAction.DATA_PROCESSING, int)
    assert isinstance(pyopenms.ProcessingAction.DEISOTOPING, int)
    assert isinstance(pyopenms.ProcessingAction.FEATURE_GROUPING, int)
    assert isinstance(pyopenms.ProcessingAction.FILTERING, int)
    assert isinstance(pyopenms.ProcessingAction.FORMAT_CONVERSION, int)
    assert isinstance(pyopenms.ProcessingAction.IDENTIFICATION_MAPPING, int)
    assert isinstance(pyopenms.ProcessingAction.NORMALIZATION, int)
    assert isinstance(pyopenms.ProcessingAction.PEAK_PICKING, int)
    assert isinstance(pyopenms.ProcessingAction.PRECURSOR_RECALCULATION, int)
    assert isinstance(pyopenms.ProcessingAction.QUANTITATION, int)
    assert isinstance(pyopenms.ProcessingAction.SIZE_OF_PROCESSINGACTION, int)
    assert isinstance(pyopenms.ProcessingAction.SMOOTHING, int)


@report
def testProduct():
    """
    @tests:
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
def testProteinHit():
    """
    @tests:
     ProteinHit.__init__
     ProteinHit.clearMetaInfo
     ProteinHit.getAccession
     ProteinHit.getCoverage
     ProteinHit.getKeys
     ProteinHit.getMetaValue
     ProteinHit.setMetaValue
     ProteinHit.getRank
     ProteinHit.__eq__
     ProteinHit.__ge__
     ProteinHit.__gt__
     ProteinHit.__le__
     ProteinHit.__lt__
     ProteinHit.__ne__
     ProteinHit.getScore
     ProteinHit.getSequence
     ProteinHit.isMetaEmpty
     ProteinHit.metaValueExists
     ProteinHit.removeMetaValue
     ProteinHit.setAccession
     ProteinHit.setCoverage
     ProteinHit.setRank
     ProteinHit.setScore
     ProteinHit.setSequence
     """
    ph = pyopenms.ProteinHit()
    assert ph == ph
    assert not ph != ph
    _testMetaInfoInterface(ph)
    ph.setAccession("A")
    ph.setCoverage(0.5)
    ph.setRank(2)
    ph.setScore(1.5)
    ph.setSequence("ABA")
    assert ph.getAccession() == ("A")
    assert ph.getCoverage() == (0.5)
    assert ph.getRank() == (2)
    assert ph.getScore() == (1.5)
    assert ph.getSequence() == ("ABA")

@report
def testProteinIdentification():
    """
    @tests:
     ProteinIdentification.DigestionEnzyme
     ProteinIdentification.PeakMassType
     ProteinIdentification.__init__
     ProteinIdentification.clearMetaInfo
     ProteinIdentification.getHits
     ProteinIdentification.getKeys
     ProteinIdentification.getMetaValue
     ProteinIdentification.insertHit
     ProteinIdentification.isMetaEmpty
     ProteinIdentification.metaValueExists
     ProteinIdentification.removeMetaValue
     ProteinIdentification.setHits
     ProteinIdentification.setMetaValue
     ProteinIdentification.__eq__
     ProteinIdentification.__ge__
     ProteinIdentification.__gt__
     ProteinIdentification.__le__
     ProteinIdentification.__lt__
     ProteinIdentification.__ne__
    """
    pi = pyopenms.ProteinIdentification()
    _testMetaInfoInterface(pi)
    assert pi == pi
    assert not pi != pi

    assert pi.getHits() == []
    ph = pyopenms.ProteinHit()
    pi.insertHit(ph)
    ph2, = pi.getHits()
    assert ph2 == ph

    pi.setHits([ph])
    ph2, = pi.getHits()
    assert ph2 == ph

    assert isinstance(pyopenms.ProteinIdentification.PeakMassType.MONOISOTOPIC, int)
    assert isinstance(pyopenms.ProteinIdentification.PeakMassType.AVERAGE, int)

    assert isinstance(pyopenms.ProteinIdentification.DigestionEnzyme.TRYPSIN,
            int)
    assert isinstance(pyopenms.ProteinIdentification.DigestionEnzyme.PEPSIN_A, int)
    assert isinstance(pyopenms.ProteinIdentification.DigestionEnzyme.PROTEASE_K,
            int)
    assert isinstance(pyopenms.ProteinIdentification.DigestionEnzyme.CHYMOTRYPSIN,
            int)
    assert isinstance(pyopenms.ProteinIdentification.DigestionEnzyme.NO_ENZYME, int)
    assert isinstance(pyopenms.ProteinIdentification.DigestionEnzyme.UNKNOWN_ENZYME,
            int)


@report
def testRichPeak():
    """
    @tests:
     RichPeak1D.__init__
     RichPeak1D.getIntensity
     RichPeak1D.getKeys
     RichPeak1D.getMZ
     RichPeak1D.__eq__
     RichPeak1D.__ge__
     RichPeak1D.__gt__
     RichPeak1D.__le__
     RichPeak1D.__lt__
     RichPeak1D.__ne__
     RichPeak1D.getMetaValue
     RichPeak1D.clearMetaInfo
     RichPeak1D.isMetaEmpty
     RichPeak1D.metaValueExists
     RichPeak1D.removeMetaValue
     RichPeak1D.setIntensity
     RichPeak1D.setMZ
     RichPeak1D.setMetaValue
     RichPeak2D.__init__
     RichPeak2D.clearUniqueId
     RichPeak2D.clearMetaInfo
     RichPeak2D.isMetaEmpty
     RichPeak2D.ensureUniqueId
     RichPeak2D.getIntensity
     RichPeak2D.getKeys
     RichPeak2D.getMZ
     RichPeak2D.getMetaValue
     RichPeak2D.getRT
     RichPeak2D.getUniqueId
     RichPeak2D.hasInvalidUniqueId
     RichPeak2D.hasValidUniqueId
     RichPeak2D.metaValueExists
     RichPeak2D.removeMetaValue
     RichPeak2D.setIntensity
     RichPeak2D.setMZ
     RichPeak2D.setMetaValue
     RichPeak2D.setUniqueId
     RichPeak2D.setRT
     RichPeak2D.__eq__
     RichPeak2D.__ge__
     RichPeak2D.__gt__
     RichPeak2D.__le__
     RichPeak2D.__lt__
     RichPeak2D.__ne__
     """
    p1 = pyopenms.RichPeak1D()
    _testMetaInfoInterface(p1)
    assert p1 == p1
    assert not p1 != p1
    p1.setMZ(12.0)
    p1.setIntensity(23.0)
    assert p1.getMZ() == (12.0)
    assert p1.getIntensity() == (23.0)

    p2 = pyopenms.RichPeak2D()
    _testMetaInfoInterface(p2)
    _testUniqueIdInterface(p2)
    assert p2 == p2
    assert not p2 != p2
    p2.setMZ(22.0)
    p2.setIntensity(23.0)
    p2.setRT(43.0)
    assert p2.getMZ() == (22.0)
    assert p2.getIntensity() == (23.0)
    assert p2.getRT() == (43.0)


@report
def testSoftware():
    """
    @tests:
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
    @tests:
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

@report
def testSpectrumSetting(s=pyopenms.SpectrumSettings()):
    """
    @tests:
     SpectrumSettings.SpectrumType
     SpectrumSettings.__init__
     SpectrumSettings.getAcquisitionInfo
     SpectrumSettings.getComment
     SpectrumSettings.getDataProcessing
     SpectrumSettings.getInstrumentSettings
     SpectrumSettings.getNativeID
     SpectrumSettings.getPeptideIdentifications
     SpectrumSettings.getPrecursors
     SpectrumSettings.getProducts
     SpectrumSettings.getSourceFile
     SpectrumSettings.getType
     SpectrumSettings.setAcquisitionInfo
     SpectrumSettings.setComment
     SpectrumSettings.setDataProcessing
     SpectrumSettings.setInstrumentSettings
     SpectrumSettings.setNativeID
     SpectrumSettings.setPeptideIdentifications
     SpectrumSettings.setPrecursors
     SpectrumSettings.setProducts
     SpectrumSettings.setSourceFile
     SpectrumSettings.setType
     SpectrumSettings.unify
    """

    assert s.getType() in [ pyopenms.SpectrumSettings.SpectrumType.UNKNOWN,
                               pyopenms.SpectrumSettings.SpectrumType.PEAKS,
                               pyopenms.SpectrumSettings.SpectrumType.RAWDATA]

    assert isinstance(s.getAcquisitionInfo(), pyopenms.AcquisitionInfo)
    assert isinstance(s.getInstrumentSettings(), pyopenms.InstrumentSettings)
    assert isinstance(s.getSourceFile(), pyopenms.SourceFile)
    assert isinstance(s.getPeptideIdentifications(), list)
    assert isinstance(s.getDataProcessing(), list)

    s.setAcquisitionInfo(s.getAcquisitionInfo())
    s.setInstrumentSettings(s.getInstrumentSettings())
    s.setSourceFile(s.getSourceFile())
    s.setPeptideIdentifications(s.getPeptideIdentifications())
    s.setDataProcessing(s.getDataProcessing())
    s.setComment(s.getComment())
    s.setPrecursors(s.getPrecursors())
    s.setProducts(s.getProducts())
    s.setType(s.getType())
    s.setNativeID(s.getNativeID())
    s.setType(s.getType())
    if isinstance(s, pyopenms.SpectrumSettings):
        s.unify(s)


@report
def testTransformationDescription():
    """
    @tests:
     TransformationDescription.__init__
     TransformationDescription.apply
     TransformationDescription.getDataPoints
     TransformationDescription.fitModel
     TransformationDescription.getModelParameters
     TransformationDescription.getModelType
     TransformationDescription.invert
    """
    td = pyopenms.TransformationDescription()
    assert td.getDataPoints() == []
    assert isinstance(td.apply(0.0), float)

    td.fitModel
    p = td.getModelParameters()
    td.getModelType()
    td.invert

@report
def testTransformationModels():
    """
    @tests:
     TransformationModelInterpolated.getDefaultParameters
     TransformationModelInterpolated.getParameters
     TransformationModelLinear.getDefaultParameters
     TransformationModelLinear.getParameters
    """
    for clz in [pyopenms.TransformationModelLinear,
                pyopenms.TransformationModelInterpolated]:
        mod = clz()
        p = pyopenms.Param()
        clz.getDefaultParameters(p)

@report
def testTransformationXMLFile():
    """
    @tests:
     TransformationXMLFile.__init__
     TransformationXMLFile.load
     TransformationXMLFile.store
    """
    fh = pyopenms.TransformationXMLFile()
    td = pyopenms.TransformationDescription()
    fh.store("test.transformationXML", td)
    fh.load("test.transformationXML", td)
    assert td.getDataPoints() == []

@report
def testType():
    """
    @tests:
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
      pyopenms.Type.CONSENSUSXML
     ,pyopenms.Type.DTA
     ,pyopenms.Type.DTA2D
     ,pyopenms.Type.EDTA
     ,pyopenms.Type.FASTA
     ,pyopenms.Type.FEATUREXML
     ,pyopenms.Type.GELML
     ,pyopenms.Type.HARDKLOER
     ,pyopenms.Type.IDXML
     ,pyopenms.Type.INI
     ,pyopenms.Type.KROENIK
     ,pyopenms.Type.MASCOTXML
     ,pyopenms.Type.MGF
     ,pyopenms.Type.MS2
     ,pyopenms.Type.MSP
     ,pyopenms.Type.MZDATA
     ,pyopenms.Type.MZIDENTML
     ,pyopenms.Type.MZML
     ,pyopenms.Type.MZXML
     ,pyopenms.Type.OMSSAXML
     ,pyopenms.Type.PEPLIST
     ,pyopenms.Type.PEPXML
     ,pyopenms.Type.PNG
     ,pyopenms.Type.PROTXML
     ,pyopenms.Type.SIZE_OF_TYPE
     ,pyopenms.Type.TOPPAS
     ,pyopenms.Type.TRAML
     ,pyopenms.Type.TRANSFORMATIONXML
     ,pyopenms.Type.TSV
     ,pyopenms.Type.UNKNOWN
     ,pyopenms.Type.XMASS]:
        assert isinstance(ti, int)

@report
def testVersion():
    """
    @tests:
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
     version.version
    """
    assert isinstance( pyopenms.VersionInfo.getVersion(), str)
    assert isinstance( pyopenms.VersionInfo.getRevision(), str)
    assert isinstance( pyopenms.VersionInfo.getTime(), str)

    vd = pyopenms.VersionDetails.create("19.2.1")
    assert vd.version_major == 19
    assert vd.version_minor == 2
    assert vd.version_patch == 1

    assert vd == vd
    assert not vd < vd
    assert not vd > vd

    assert  isinstance(pyopenms.version.version, str)


@report
def testPILISCrossValidation():
    """
    @tests:
     PILISCrossValidation.__init__
    """
    inst = pyopenms.PILISCrossValidation()

    assert inst.apply is not None
    assert inst.setOption is not None

@report
def testPILIS_Peptide():
    """
    @tests:
     PILIS_Peptide.__init__
    """
    inst = pyopenms.PILIS_Peptide()

    assert inst.sequence is not None
    assert inst.charge is not None
    assert inst.spec is not None
    assert inst.hits is not None
    
@report
def testPILIS_Option():
    """
    @tests:
     PILIS_Option.__init__
    """
    inst = pyopenms.PILIS_Option()

    assert inst.type is not None
    assert inst.int_min is not None
    assert inst.int_max is not None
    assert inst.int_stepsize is not None
    assert inst.dbl_min is not None
    assert inst.dbl_max is not None
    assert inst.dbl_stepsize is not None

@report
def testPILIS_Option_Type():
    """
    @tests:
     PILIS_Option_Type.__init__
    """
    inst = pyopenms.PILIS_Option.PILIS_Option_Type()

    assert inst.INT is not None
    assert inst.DOUBLE is not None
    assert inst.BOOL is not None
    assert inst.STRINGLIST is not None

@report
def testInspectInfile():
    """
    @tests:
     InspectInfile.__init__
    """
    inst = pyopenms.InspectInfile()

    assert inst.getModifications is not None
    mods = inst.getModifications()
    assert len(mods) == 0


@report
def testIsotopeMarker():
    """
    @tests:
     IsotopeMarker.__init__
    """
    inst = pyopenms.IsotopeMarker()
    ptr = inst.create()

    assert ptr.apply is not None

    res = {}
    spec = pyopenms.MSSpectrum()
    ptr.apply(res, spec)

@report
def testAttachment():
    """
    @tests:
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
    inst.colTypes = [ "test", "test2"]
    inst.tableRows = [ ["test", "test2"], ["otherTest"] ]

    assert inst.tableRows[1][0] == "otherTest"

@report
def testOptimizePeakDeconvolution():
    """
    @tests:
     OptimizePeakDeconvolution.__init__
    """
    inst = pyopenms.OptimizePeakDeconvolution()
    assert inst.getParameters

    assert inst.getPenalties is not None
    assert inst.setPenalties is not None
    assert inst.getCharge is not None
    assert inst.setCharge is not None
    assert inst.optimize is not None


    inst = pyopenms.PenaltyFactorsIntensity()
    assert inst.height is not None

    inst = pyopenms.OptimizePeakDeconvolution_Data()
    assert inst.peaks is not None
    assert inst.peaks is not None
    assert inst.signal is not None
    assert inst.penalties is not None
    assert inst.charge is not None


@report
def testIndexedMzMLDecoder():
    decoder = pyopenms.IndexedMzMLDecoder()
    pos = decoder.findIndexListOffset("abcde", 100)
    assert isinstance(pos, pyopenms.streampos)
    assert long(pos) == -1   # not found


def test_streampos():
    p = pyopenms.streampos()
    assert isinstance(int(p), int)

def test_MapConversion():

    feature = pyopenms.Feature()
    feature.setRT(99)

    cmap = pyopenms.ConsensusMap()
    fmap = pyopenms.FeatureMap()
    fmap.push_back(feature)
    pyopenms.MapConversion().convert(0, fmap, cmap, 1)

    assert(cmap.size() == 1)
    assert(cmap[0].getRT() == 99.0)

    fmap = pyopenms.FeatureMap()
    pyopenms.MapConversion().convert(cmap, True, fmap)

    assert(fmap.size() == 1)
    assert(fmap[0].getRT() == 99.0)

    exp = pyopenms.MSExperiment()
    sp = pyopenms.MSSpectrum()
    peak = pyopenms.Peak1D()
    peak.setIntensity(10)
    peak.setMZ(20)
    sp.push_back(peak)
    exp.addSpectrum(sp)
    exp.addSpectrum(sp)

    cmap = pyopenms.ConsensusMap()
    pyopenms.MapConversion().convert(0, exp, cmap, 2)

    assert(cmap.size() == 2)
    assert(cmap[0].getIntensity() == 10.0)
    assert(cmap[0].getMZ() == 20.0)

