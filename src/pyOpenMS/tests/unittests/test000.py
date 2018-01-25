#!/usr/bin/env python
# -*- coding: utf-8  -*-
from __future__ import print_function

import pyopenms
import copy

from pyopenms import String as s

print(b"IMPORTED b", pyopenms.__file__)

try:
    long
except NameError:
    long = int

from functools import wraps


def report(f):
    @wraps(f)
    def wrapper(*a, **kw):
        print(b"run b", f.__name__)
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

    what.setMetaValue(b"key", 42)
    what.setMetaValue(b"key2", 42)

    keys = []
    what.getKeys(keys)
    assert len(keys) and all(isinstance(k, bytes) for k in keys)
    assert what.getMetaValue(keys[0]) == 42
    keys = []
    what.getKeysAsIntegers(keys)
    assert len(keys) and all(isinstance(k, (long, int)) for k in keys)

    assert what.metaValueExists(b"key")
    what.removeMetaValue(b"key")

    what.setMetaValue(1024, 42)
    assert what.getMetaValue(1024) == 42

    keys = []
    what.getKeys(keys)
    assert what.getMetaValue(keys[0]) == 42
    keys = []
    what.getKeysAsIntegers(keys)
    assert len(keys) and all(isinstance(k, (long, int)) for k in keys)

    what.clearMetaInfo()
    keys = []
    what.getKeys(keys)
    assert len(keys) == 0
    what.getKeysAsIntegers(keys)
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
    ff.startProgress(0, 3, b"label")
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

    rich_spec = pyopenms.MSSpectrum()
    p = pyopenms.Peak1D()
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
     AASequence.getCTerminalModificationName
     AASequence.getNTerminalModificationName
     AASequence.setCTerminalModification
     AASequence.setModification
     AASequence.setNTerminalModification
     AASequence.toString
     AASequence.toUnmodifiedString
    """
    aas = pyopenms.AASequence()

    aas + aas
    aas += aas

    aas.__doc__
    aas = pyopenms.AASequence.fromString(b"DFPIANGER", True)
    assert aas.getCTerminalModificationName() == b""
    assert aas.getNTerminalModificationName() == b""
    aas.setCTerminalModification(b"")
    aas.setNTerminalModification(b"")
    assert aas.toString() == b"DFPIANGER"
    assert aas.toUnmodifiedString() == b"DFPIANGER"

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
    ins.setName(b"Carbon")
    ins.getName()
    ins.setSymbol(b"C")
    ins.getSymbol()

    e = pyopenms.Element()
    e.setSymbol("blah")
    e.setSymbol(b"blah")
    e.setSymbol(u"blah")
    e.setSymbol(str("blah"))
    oms_string = s(b"blub")
    e.setSymbol(oms_string)
    assert oms_string
    assert oms_string.toString() == "blub"

    evil = u"blüb"
    evil8 = evil.encode("utf8")
    evil1 = evil.encode("latin1")


    e.setSymbol(evil.encode("utf8"))
    assert e.getSymbol().decode("utf8") == u"blüb"
    e.setSymbol(evil.encode("latin1"))
    assert e.getSymbol().decode("latin1") == u"blüb"

    # If we get the raw symbols, we get bytes (which we would need to decode first)
    e.setSymbol(evil8.decode("utf8"))
    assert e.getSymbol() == b'bl\xc3\xbcb'
    assert e.getSymbol() == u"blüb".encode("utf8")
    # OpenMS strings, however, understand the decoding
    assert s(e.getSymbol()) == s(u"blüb")
    assert s(e.getSymbol()).toString() == u"blüb"

    # What if you use the wrong decoding ?
    e.setSymbol(evil1)
    #print(e.getSymbol().decode("latin1"))
    assert e.getSymbol().decode("latin1") == u"blüb"
    e.setSymbol(evil8)
    #print(e.getSymbol().decode("utf8"))
    assert e.getSymbol().decode("latin1") == u"blüb".encode("utf8").decode("latin1")

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
    pyopenms.Residue.ResidueType.CIon
    pyopenms.Residue.ResidueType.XIon
    pyopenms.Residue.ResidueType.YIon
    pyopenms.Residue.ResidueType.ZIon
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
    # ins.getNumberOf(b"test")
    ins.getNumberOfAtoms()
    ins.setCharge(2)
    ins.getCharge()
    ins.toString()
    ins.isEmpty()
    ins.isCharged()
    ins.hasElement( pyopenms.Element() )

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

    f.setId(b"test_id")
    assert f.getId() == b"test_id"

    f.setId("test_id")
    assert f.getId() == b"test_id"

    f.setCharge(5)
    assert f.getCharge() == 5

    f.setCalculatedMassToCharge(5.0)
    assert f.getCalculatedMassToCharge() == 5.0

    f.setExperimentalMassToCharge(5.0)
    assert f.getExperimentalMassToCharge() == 5.0

    f.setName(b"test")
    assert f.getName() == b"test"

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
    hit.setName(b"test1")
    f.addHit(hit)
    hit = pyopenms.IdentificationHit()
    hit.setName(b"test2")
    f.addHit(hit)
    all_hits = f.getHits()
    assert len(all_hits) == 2
    assert b"test1" in [h.getName() for h in all_hits]
    assert b"test2" in [h.getName() for h in all_hits]

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
    fixed = [b"Carbamidomethyl"]
    variable = [b"Oxidation"]
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
    ai.setMethodOfCombination(b"ABC")
    assert ai.getMethodOfCombination() == b"ABC"

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
    f_ = copy.copy(f)
    assert f_ == f
    f_ = copy.deepcopy(f)
    assert f_ == f
    f_ = pyopenms.ConsensusFeature(f)
    assert f_ == f

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

@report
def testXTandemInfile():
    """
    """
    f = pyopenms.XTandemInfile()

    f.setFragmentMassTolerance is not None
    f.getFragmentMassTolerance is not None

    f.setPrecursorMassTolerancePlus is not None
    f.getPrecursorMassTolerancePlus is not None
    f.setPrecursorMassToleranceMinus is not None
    f.getPrecursorMassToleranceMinus is not None

    f.setPrecursorErrorType is not None
    f.getPrecursorErrorType is not None

    f.setFragmentMassErrorUnit is not None
    f.getFragmentMassErrorUnit is not None
    f.setPrecursorMassErrorUnit is not None
    f.getPrecursorMassErrorUnit is not None

    f.setNumberOfThreads is not None
    f.getNumberOfThreads is not None

    f.setModifications is not None
    f.getModifications is not None

    f.setOutputFilename is not None
    f.getOutputFilename is not None
    f.setInputFilename is not None
    f.getInputFilename is not None
    f.setTaxonomyFilename is not None
    f.getTaxonomyFilename is not None
    f.setDefaultParametersFilename is not None
    f.getDefaultParametersFilename is not None


    f.setTaxon(b"testTaxon")
    assert f.getTaxon() == b"testTaxon"

    assert f.setMaxPrecursorCharge is not None
    assert f.getMaxPrecursorCharge is not None

    assert f.setNumberOfMissedCleavages is not None
    assert f.getNumberOfMissedCleavages is not None

    assert f.setMaxValidEValue is not None
    assert f.getMaxValidEValue is not None

    assert f.setSemiCleavage is not None
    assert f.setAllowIsotopeError is not None
    assert f.write is not None
    assert f.setCleavageSite is not None
    assert f.getCleavageSite is not None

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

    assert isinstance(dp.getCompletionTime().getDate(), bytes)
    assert isinstance(dp.getCompletionTime().getTime(), bytes)
    dp.clearMetaInfo()
    k = []
    dp.getKeys(k)
    assert k == []
    dp.getMetaValue
    ac = dp.getProcessingActions()
    assert ac == set(())
    ac = set([ pyopenms.ProcessingAction.PEAK_PICKING, pyopenms.ProcessingAction.BASELINE_REDUCTION])
    dp.setProcessingActions(ac)
    assert len(dp.getProcessingActions() ) == 2
    assert isinstance(dp.getSoftware().getName(), bytes)
    assert isinstance(dp.getSoftware().getVersion(), bytes)
    dp.isMetaEmpty()
    dp.metaValueExists
    dp.removeMetaValue
    dp.setCompletionTime(pyopenms.DateTime.now())
    s = dp.getSoftware()
    s.setName(b"pyopenms")
    dp.setSoftware(s)

    assert dp.getSoftware().getName() == b"pyopenms"


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

    a = pyopenms.DataValue(b"1")
    assert not a.isEmpty()
    assert a.toString() == b"1"
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

    assert pyopenms.MSSpectrum().getMetaValue(b"nonexisingkey") is None

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
    assert isinstance( d.getDate(), bytes)
    assert isinstance( d.getTime(), bytes)
    d = pyopenms.DateTime.now()
    assert isinstance( d.getDate(), bytes)
    assert isinstance( d.getTime(), bytes)

    d.clear()
    d.set(b"01.01.2001 11:11:11")
    assert d.get() == b"2001-01-01 11:11:11"

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
        if len(k.split(b":")) < 2: continue
        f = k.split(b":")[0]
        p.setSectionDescription(f, k)
        assert p.getSectionDescription(f) == k

        assert p.get(k) is not None

    assert len(p.values()) == len([p[k] for k in p.keys()])
    assert sorted(p.items()) == sorted((k, p[k]) for k in p.keys())

    assert not p.exists(b"asdflkj01231321321v")
    p.addTag(k, b"a")
    p.addTags(k, [b"b", b"c"])
    assert sorted(p.getTags(k)) == [b"a", b"b", b"c"]
    p.clearTags(k)
    assert p.getTags(k) == []

    pn = pyopenms.Param()
    pn.insert(b"master:", p)
    assert pn.exists(b"master:"+k)

    p1 = pn.copy(b"master:", True)
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
    ph.store(b"test.ini", p)
    p1 = pyopenms.Param()
    ph.load(b"test.ini", p1)
    assert p == p1

    e1 = p1.getEntry(k)
    for f in [b"name", b"description", b"value", b"tags", b"valid_strings",
              b"min_float", b"max_float", b"min_int", b"max_int"]:
        assert getattr(e1, f.decode()) is not None

    assert e1 == e1

    assert p1.get(b"abcde", 7) == 7




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

    assert ff.getName() == b"FeatureFinderAlgorithm"
    assert pyopenms.FeatureFinderAlgorithmPicked.getProductName() == b"centroided"

    ff.setParameters(pyopenms.Param())

    ff.setName(b"test")
    assert ff.getName() == b"test"

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

    assert ff.getName() == b"FeatureFinderAlgorithm"
    assert pyopenms.FeatureFinderAlgorithmSH.getProductName() == b"superhirn"

    ff.setParameters(pyopenms.Param())

    ff.setName(b"test")
    assert ff.getName() == b"test"

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

    assert ff.getName() == b"FeatureFinderAlgorithm"
    assert pyopenms.FeatureFinderAlgorithmIsotopeWavelet.getProductName() == b"isotope_wavelet"

    ff.setParameters(pyopenms.Param())

    ff.setName(b"test")
    assert ff.getName() == b"test"

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

# TODO: re-enable as soon as ConsensusIDAlgorithm classes are wrapped
# @report
# def testConsensusID():
#     """
#     @tests:
#      ConsensusID.__init__
#     """
#     ff = pyopenms.ConsensusID()
#     p = ff.getDefaults()
#     _testParam(p)

#     assert pyopenms.ConsensusID().apply is not None

@report
def testFalseDiscoveryRate():
    """
    @tests:
     FalseDiscoveryRate.__init__
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
    spectrum = pyopenms.MSSpectrum()

    ff.compute(hit, spectrum)
    # ff.computeCumulativeScore(1,1,0.5)

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

    assert pyopenms.InternalCalibration().fillCalibrants is not None
    assert pyopenms.InternalCalibration().getCalibrationPoints is not None
    assert pyopenms.InternalCalibration().calibrate is not None

@report
def testTransitionTSVFile():
    """
    @tests:
     TransitionTSVFile.__init__
     TransitionTSVFile.calibrateMapGlobally
     TransitionTSVFile.calibrateMapSpectrumwise
    """
    ff = pyopenms.TransitionTSVFile()

    assert pyopenms.TransitionTSVFile().convertTargetedExperimentToTSV is not None
    assert pyopenms.TransitionTSVFile().convertTSVToTargetedExperiment is not None
    assert pyopenms.TransitionTSVFile().validateTargetedExperiment is not None

@report
def testProteaseDigestion():
    """
    @tests:
     ProteaseDigestion.__init__
     ProteaseDigestion.getMissedCleavages()
     ProteaseDigestion.setMissedCleavages()
     ProteaseDigestion.digest()
     ProteaseDigestion.peptideCount()
    """
    # removed due to name clashes
    # ProteaseDigestion.getEnzyme()
    # ProteaseDigestion.setEnzyme()
    # ProteaseDigestion.getEnzymeByName()

    ff = pyopenms.ProteaseDigestion()
    #enz = pyopenms.ProteaseDigestion().Enzyme()

    assert pyopenms.ProteaseDigestion().getMissedCleavages is not None
    assert pyopenms.ProteaseDigestion().setMissedCleavages is not None
    #assert pyopenms.ProteaseDigestion().getEnzyme is not None
    #assert pyopenms.ProteaseDigestion().setEnzyme is not None
    #assert pyopenms.ProteaseDigestion().getEnzymeByName is not None

    assert pyopenms.ProteaseDigestion().digest is not None
    assert pyopenms.ProteaseDigestion().peptideCount is not None

    ff.setMissedCleavages(5)
    assert ff.getMissedCleavages() == 5

    #ff.setEnzyme(enz.TRYPSIN)
    #assert ff.getEnzyme() == enz.TRYPSIN

@report
def testEnzymaticDigestionLogModel():
    ff = pyopenms.EnzymaticDigestionLogModel()
    assert pyopenms.EnzymaticDigestionLogModel().getLogThreshold is not None
    assert pyopenms.EnzymaticDigestionLogModel().setLogThreshold is not None
    assert pyopenms.EnzymaticDigestionLogModel().digest is not None
    assert pyopenms.EnzymaticDigestionLogModel().peptideCount is not None
    ff.setLogThreshold(0.25)
    assert ff.getLogThreshold() == 0.25

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
    fm_ = copy.copy(fm)
    assert fm_ == fm
    fm_ = copy.deepcopy(fm)
    assert fm_ == fm
    fm_ = pyopenms.FeatureMap(fm)
    assert fm_ == fm

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
    assert fm + fm2 != fm


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
    fh.store(b"test.featureXML", fm)
    fh.load(b"test.featureXML", fm)

    fh = pyopenms.FileHandler()
    fh.loadFeatures(b"test.featureXML", fm)

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
    assert isinstance(fd.filename, bytes)
    assert isinstance(fd.label, bytes)
    assert isinstance(fd.size, int)
    assert isinstance(fd.unique_id, (long, int, bytes))

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
    fh.storeExperiment(b"test1.mzML", mse)
    fh.loadExperiment(b"test1.mzML", mse)
    fh.storeExperiment(b"test1.mzXML", mse)
    fh.loadExperiment(b"test1.mzXML", mse)
    fh.storeExperiment(b"test1.mzData", mse)
    fh.loadExperiment(b"test1.mzData", mse)

    

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
    idm.setName(b"x")
    assert idm.getName() == b"x"
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
    ins.setFirstName(b"test")
    ins.getLastName()
    ins.setLastName(b"test")
    ins.setName(b"Testy Test")
    ins.getInstitution()
    ins.setInstitution(b"test")
    ins.getEmail()
    ins.setEmail(b"test")
    ins.getURL()
    ins.setURL(b"test")
    ins.getAddress()
    ins.setAddress(b"test")
    ins.getContactInfo()
    ins.setContactInfo(b"test")

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

    ins.setIdentifier(b"test")
    ins.getIdentifier()
    # ins.setLoadedFilePath(b"Test")
    ins.getLoadedFilePath()
    # ins.setLoadedFileType(b"test")
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

    ins.addEluent(b"test")
    ins.clearEluents()
    assert len(ins.getEluents() ) == 0
    ins.addEluent(b"test")
    assert len(ins.getEluents() ) == 1

    ins.clearTimepoints()
    ins.addTimepoint(5)
    ins.getTimepoints()

    ins.setPercentage(b"test", 5, 20)
    ins.getPercentage(b"test", 5)
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

    ins.setInstrument(b"test")
    ins.getInstrument()
    ins.setColumn(b"test")
    ins.getColumn()
    ins.setTemperature(6)
    ins.getTemperature()
    ins.setPressure(6)
    ins.getPressure()
    ins.setFlux(8)
    ins.getFlux()
    ins.setComment(b"test")
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

    ins.setName(b"test")
    ins.getName()
    ins.setVendor(b"test")
    ins.getVendor()
    ins.setModel(b"test")
    ins.getModel()
    ins.setCustomizations(b"test")
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

    ins.setName(b"test")
    ins.getName()
    ins.setOrganism(b"test")
    ins.getOrganism()
    ins.setNumber(b"test")
    ins.getNumber()
    ins.setComment(b"test")
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
    mse_ = copy.copy(mse)
    assert mse_ == mse
    mse_ = copy.deepcopy(mse)
    assert mse_ == mse
    mse_ = pyopenms.MSExperiment(mse)
    assert mse_ == mse

    _testMetaInfoInterface(mse)
    mse.updateRanges()
    mse.sortSpectra(True)
    assert isinstance(mse.getMaxRT(), float)
    assert isinstance(mse.getMinRT(), float)
    assert isinstance(mse.getMaxMZ(), float)
    assert isinstance(mse.getMinMZ(), float)
    assert isinstance(mse.getLoadedFilePath(), bytes)
    assert isinstance(mse.getMinInt(), float)
    assert isinstance(mse.getMaxInt(), float)

    assert isinstance(mse.getMin()[0], float)
    assert isinstance(mse.getMin()[1], float)
    assert isinstance(mse.getMax()[0], float)
    assert isinstance(mse.getMax()[1], float)
    mse.setLoadedFilePath(b"")
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
    spec_ = copy.copy(spec)
    assert spec_ == spec
    spec_ = copy.deepcopy(spec)
    assert spec_ == spec
    spec_ = pyopenms.MSSpectrum(spec)
    assert spec_ == spec

    _testMetaInfoInterface(spec)

    testSpectrumSetting(spec)

    spec.setRT(3.0)
    assert spec.getRT() == 3.0
    spec.setMSLevel(2)
    assert spec.getMSLevel() == 2
    spec.setName(b"spec")
    assert spec.getName() == b"spec"

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

    mz, ii = spec.get_peaks()
    assert len(mz) == len(ii)
    assert len(mz) == 1

    spec.set_peaks((mz, ii))
    mz0, ii0 = spec.get_peaks()
    assert mz0 == mz
    assert ii0 == ii

    assert int(spec.isSorted()) in  (0,1)

    # get data arrays
    assert len(spec.getStringDataArrays()) == 0
    string_da = [ pyopenms.StringDataArray() ]
    string_da[0].push_back("hello")
    string_da[0].push_back("world")
    string_da.append( pyopenms.StringDataArray() )
    string_da[1].push_back("other")
    spec.setStringDataArrays( string_da )
    assert len(spec.getStringDataArrays()) == 2
    assert spec.getStringDataArrays()[0][0] == b"hello"
    assert spec.getStringDataArrays()[1][0] == b"other"

    assert len(spec.getIntegerDataArrays()) == 0
    # int_da = [ [5, 6], [8] ]
    int_da = [ pyopenms.IntegerDataArray() ]
    int_da[0].push_back(5)
    int_da[0].push_back(6)
    int_da.append( pyopenms.IntegerDataArray() )
    int_da[1].push_back(8)
    spec.setIntegerDataArrays( int_da )
    assert len(spec.getIntegerDataArrays()) == 2
    assert spec.getIntegerDataArrays()[0][0] == 5
    assert spec.getIntegerDataArrays()[1][0] == 8

    assert len(spec.getFloatDataArrays()) == 0
    # int_da = [ [5, 6], [8] ]
    int_da = [ pyopenms.FloatDataArray() ]
    int_da[0].push_back(5.0)
    int_da[0].push_back(6.0)
    int_da.append( pyopenms.FloatDataArray() )
    int_da[1].push_back(8.0)
    spec.setFloatDataArrays( int_da )
    assert len(spec.getFloatDataArrays()) == 2.0
    assert spec.getFloatDataArrays()[0][0] == 5.0
    assert spec.getIntegerDataArrays()[1][0] == 8

@report
def testStringDataArray():
    """
    @tests:
     """
    da = pyopenms.StringDataArray()
    assert da.size() == 0
    da.push_back("hello")
    da.push_back("world")
    assert da.size() == 2
    assert da[0] == b"hello"
    assert da[1] == b"world"
    da[1] = "hello world"
    assert da[1] == b"hello world", da[1]
    da.clear()
    assert da.size() == 0
    da.push_back("hello")
    assert da.size() == 1
    da.resize(3)
    da[0] = "hello"
    da[1] = ""
    da[2] = "world"
    assert da.size() == 3

@report
def testIntegerDataArray():
    """
    @tests:
     """
    da = pyopenms.IntegerDataArray()
    assert da.size() == 0
    da.push_back(1)
    da.push_back(4)
    assert da.size() == 2
    assert da[0] == 1
    assert da[1] == 4
    da[1] = 7
    assert da[1] == 7
    da.clear()
    assert da.size() == 0
    da.push_back(1)
    assert da.size() == 1
    da.resize(3)
    da[0] = 1
    da[1] = 2
    da[2] = 3
    assert da.size() == 3

@report
def testFloatDataArray():
    """
    @tests:
     """
    da = pyopenms.FloatDataArray()
    assert da.size() == 0
    da.push_back(1.0)
    da.push_back(4.0)
    assert da.size() == 2
    assert da[0] == 1.0
    assert da[1] == 4.0
    da[1] = 7.0
    assert da[1] == 7.0
    da.clear()
    assert da.size() == 0
    da.push_back(1.0)
    assert da.size() == 1
    da.resize(3)
    da[0] = 1.0
    da[1] = 2.0
    da[2] = 3.0
    assert da.size() == 3

@report
def testMSChromatogram():
    """
    @tests:
     MSChromatogram.__init__
     MSChromatogram.__copy__
     """
    chrom = pyopenms.MSChromatogram()
    chrom_ = copy.copy(chrom)
    assert chrom_ == chrom
    chrom_ = copy.deepcopy(chrom)
    assert chrom_ == chrom
    chrom_ = pyopenms.MSChromatogram(chrom)
    assert chrom_ == chrom

    _testMetaInfoInterface(chrom)

    chrom.setName(b"chrom")
    assert chrom.getName() == b"chrom"

    p = pyopenms.ChromatogramPeak()
    p.setRT(1000.0)
    p.setIntensity(200.0)

    chrom.push_back(p)
    assert chrom.size() == 1
    assert chrom[0] == p

    chrom.updateRanges()
    assert isinstance(chrom.findNearest(0.0), int)

    assert isinstance(chrom.getMin()[0], float)
    assert isinstance(chrom.getMax()[0], float)
    assert isinstance(chrom.getMinInt(), float)
    assert isinstance(chrom.getMaxInt(), float)

    assert chrom == chrom
    assert not chrom != chrom

    mz, ii = chrom.get_peaks()
    assert len(mz) == len(ii)
    assert len(mz) == 1

    chrom.set_peaks((mz, ii))
    mz0, ii0 = chrom.get_peaks()
    assert mz0 == mz
    assert ii0 == ii

    assert int(chrom.isSorted()) in  (0,1)

@report
def testMRMFeature():
    """
    @tests:
      MRMFeature.__init__
      MRMFeature.addScore
      MRMFeature.getScore
     """
    mrmfeature = pyopenms.MRMFeature()

    mrmfeature.addScore(b"testscore", 6)
    assert mrmfeature.getScore(b"testscore") == 6.0
    mrmfeature.addScore(b"testscore", 7)
    assert mrmfeature.getScore(b"testscore") == 7.0

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

    assert pyopenms.MRMDecoy().generateDecoys is not None

@report
def testMRMTransitionGroup():
    """
    @tests:
     """
    mrmgroup = pyopenms.MRMTransitionGroupCP()
    assert mrmgroup is not None

    mrmgroup.setTransitionGroupID(b"this_id")
    assert mrmgroup.getTransitionGroupID() == b"this_id"

    assert len(mrmgroup.getTransitions()) == 0
    mrmgroup.addTransition(pyopenms.ReactionMonitoringTransition(), b"tr1")
    assert len(mrmgroup.getTransitions()) == 1

@report
def testReactionMonitoringTransition():
    """
    @tests:
     """
    tr = pyopenms.ReactionMonitoringTransition()

@report
def testTargetedExperiment():
    """
    @tests: TargetedExperiment
     """
    m = pyopenms.TargetedExperiment()
    m_ = copy.copy(m)
    assert m_ == m
    m_ = copy.deepcopy(m)
    assert m_ == m
    m_ = pyopenms.TargetedExperiment(m)
    assert m_ == m

    m.clear(True)
    m.setCVs(m.getCVs())

    targeted = m

    targeted.setCVs(targeted.getCVs())
    targeted.setTargetCVTerms(targeted.getTargetCVTerms())
    targeted.setPeptides(targeted.getPeptides())
    targeted.setProteins(targeted.getProteins())
    targeted.setTransitions(targeted.getTransitions())

    assert m == m
    assert not m != m


@report
def testTargetedExperimentHelper():
    """
    @tests: TargetedExperimentHelper
     """
    rtu = pyopenms.RetentionTime.RTUnit()
    rtu = pyopenms.RetentionTime.RTUnit.SECOND
    rtu = pyopenms.RetentionTime.RTUnit.MINUTE
    rtt = pyopenms.RetentionTime.RTType()
    rtt = pyopenms.RetentionTime.RTType.LOCAL
    rtt = pyopenms.RetentionTime.RTType.NORMALIZED
    rtt = pyopenms.RetentionTime.RTType.IRT

    rt = pyopenms.RetentionTime()
    assert rt.software_ref is not None
    assert not rt.isRTset()
    rt.setRT(5.0)
    rt.retention_time_unit = pyopenms.RetentionTime.RTUnit.SECOND
    rt.retention_time_type = pyopenms.RetentionTime.RTType.NORMALIZED
    assert rt.isRTset()
    assert rt.getRT() == 5.0

    p = pyopenms.Peptide()
    assert p.rts is not None
    assert p.id is not None
    assert p.protein_refs is not None
    assert p.evidence is not None
    assert p.sequence is not None
    assert p.mods is not None

    assert not p.hasCharge()
    p.setChargeState(5)
    assert p.hasCharge()
    assert p.getChargeState() == 5

    assert not p.hasRetentionTime()
    p.rts = [rt]
    assert p.hasRetentionTime()
    assert p.getRetentionTime() == 5.0
    assert p.getRetentionTimeUnit() == pyopenms.RetentionTime.RTUnit.SECOND
    assert p.getRetentionTimeType() == pyopenms.RetentionTime.RTType.NORMALIZED

    c = pyopenms.Compound()
    assert c.rts is not None
    assert c.id is not None
    assert c.molecular_formula is not None
    assert c.smiles_string is not None
    assert c.theoretical_mass is not None

    assert not c.hasCharge()
    c.setChargeState(5)
    assert c.hasCharge()
    assert c.getChargeState() == 5

    assert not c.hasRetentionTime()
    c.rts = [rt]
    assert c.hasRetentionTime()
    assert c.getRetentionTime() == 5.0
    assert c.getRetentionTimeUnit() == pyopenms.RetentionTime.RTUnit.SECOND
    assert c.getRetentionTimeType() == pyopenms.RetentionTime.RTType.NORMALIZED

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

     MapAlignmentTransformer.transformRetentionTimes
     """
    ma = pyopenms.MapAlignmentAlgorithmPoseClustering()
    assert isinstance(ma.getDefaults(), pyopenms.Param)
    assert isinstance(ma.getParameters(), pyopenms.Param)
    assert isinstance(ma.getName(), bytes)

    ma.setName(ma.getName())

    ma.getDefaults()
    ma.getParameters()

    ma.setParameters(ma.getDefaults())

    ma.setReference
    ma.align

    pyopenms.MapAlignmentTransformer.transformRetentionTimes

@report
def testMapAlignmentIdentification():

    """
    @tests:
     MapAlignmentAlgorithmIdentification.__init__
     """
    ma = pyopenms.MapAlignmentAlgorithmIdentification()

    assert pyopenms.MapAlignmentAlgorithmIdentification().align is not None
    assert pyopenms.MapAlignmentAlgorithmIdentification().setReference is not None

@report
def testMapAlignmentTransformer():

    """
    @tests:
     MapAlignmentTransformer.__init__
     """
    ma = pyopenms.MapAlignmentTransformer()

    assert pyopenms.MapAlignmentTransformer().transformRetentionTimes is not None

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
    fh.store(b"test.mzData", mse)
    fh.load(b"test.mzData", mse)

    fh.setOptions(fh.getOptions())

    fh = pyopenms.MzMLFile()
    _testProgressLogger(fh)
    fh.store(b"test.mzML", mse)
    fh.load(b"test.mzML", mse)
    fh.setOptions(fh.getOptions())

    fh = pyopenms.MzXMLFile()
    _testProgressLogger(fh)
    fh.store(b"test.mzXML", mse)
    fh.load(b"test.mzXML", mse)
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
    fh.store(b"test.ini", p)
    fh.load(b"test.ini", p)



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
def testMRMMapping():
    """
    @tests:
     MRMMapping.__init__
     MRMMapping.map
    """

    p = pyopenms.MRMMapping()
    assert p.mapExperiment is not None
    e = pyopenms.MSExperiment()
    c = pyopenms.MSChromatogram()
    e.addChromatogram(c)
    assert e.getNrChromatograms() == 1

    o = pyopenms.MSExperiment()
    t = pyopenms.TargetedExperiment()
    p.mapExperiment(e, t, o)
    assert o.getNrChromatograms() == 0 # not so easy to test

@report
def testPeakPickerMRM():
    """
    @tests:
     PeakPickerMRM.__init__
     PeakPickerMRM.pickChromatogram
    """

    p = pyopenms.PeakPickerMRM()
    assert p.pickChromatogram is not None

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

    p = pyopenms.PeakPickerHiRes()
    assert p.pick is not None
    assert p.pickExperiment is not None

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

    ph = pyopenms.PeptideHit(1.0, 1, 0, pyopenms.AASequence.fromString(b"A", True))
    _testMetaInfoInterface(ph)

    assert len(ph.getPeptideEvidences()) == 0
    assert ph.getPeptideEvidences() == []

    pe = pyopenms.PeptideEvidence()
    pe.setProteinAccession(b'B_id')

    ph.addPeptideEvidence(pe)
    assert len(ph.getPeptideEvidences()) == 1
    assert ph.getPeptideEvidences()[0].getProteinAccession() == b'B_id'

    ph.setPeptideEvidences([pe,pe])
    assert len(ph.getPeptideEvidences()) == 2
    assert ph.getPeptideEvidences()[0].getProteinAccession() == b'B_id'

    assert ph.getScore() == 1.0
    assert ph.getRank() == 1
    assert ph.getSequence().toString() == b"A"

    ph.setScore(2.0)
    assert ph.getScore() == 2.0
    ph.setRank(30)
    assert ph.getRank() == 30
    ph.setSequence(pyopenms.AASequence.fromString(b"AAA", True))
    assert ph.getSequence().toString() == b"AAA"

    assert ph == ph
    assert not ph != ph

@report
def testPeptideEvidence():
    """
    @tests:
     PeptideEvidence.__init__
    """
    pe = pyopenms.PeptideEvidence()
    assert pe == pe
    assert not pe != pe

    pe.setProteinAccession(b'B_id')
    assert pe.getProteinAccession() == b"B_id"

    pe.setAABefore(b'B')
    assert pe.getAABefore() == 'B'
    pe.setAAAfter(b'C')
    assert pe.getAAAfter() == 'C'

    pe.setStart(5)
    assert pe.getStart() == 5
    pe.setEnd(9)
    assert pe.getEnd() == 9

    assert pe == pe
    assert not pe != pe


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

    pe = pyopenms.PeptideEvidence()
    pe.setProteinAccession(b'B_id')

    ph = pyopenms.PeptideHit(1.0, 1, 0, pyopenms.AASequence.fromString(b"A", True))
    ph.addPeptideEvidence(pe)
    pi.insertHit(ph)
    phx, = pi.getHits()
    assert phx == ph

    pi.setHits([ph])
    phx, = pi.getHits()
    assert phx == ph

    rv = set([])
    peptide_hits = pi.getReferencingHits(pi.getHits(), rv)
    assert rv == set([])
    # assert len(peptide_hits) == 1

    assert isinstance(pi.getSignificanceThreshold(), float)
    assert isinstance(pi.getScoreType(), bytes)
    pi.setScoreType(b"A")
    assert isinstance(pi.isHigherScoreBetter(), int)
    assert isinstance(pi.getIdentifier(), bytes)
    pi.setIdentifier(b"id")
    pi.assignRanks()
    pi.sort()
    assert not pi.empty()

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
    ph.setAccession(b"A")
    ph.setCoverage(0.5)
    ph.setRank(2)
    ph.setScore(1.5)
    ph.setSequence(b"ABA")
    assert ph.getAccession() == (b"A")
    assert ph.getCoverage() == (0.5)
    assert ph.getRank() == (2)
    assert ph.getScore() == (1.5)
    assert ph.getSequence() == (b"ABA")

@report
def testProteinIdentification():
    """
    @tests:
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
    sw.setName(b"name")
    sw.setVersion(b"1.0.0")
    assert sw.getName() == b"name"
    assert sw.getVersion() == b"1.0.0"



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
    sf.setNameOfFile(b"file.txt")
    assert sf.getNameOfFile() == b"file.txt"
    sf.setPathToFile(b"file.txt")
    assert sf.getPathToFile() == b"file.txt"
    sf.setFileType(b".txt")
    assert sf.getFileType() == b".txt"
    sf.setChecksum(b"abcde000", pyopenms.ChecksumType.UNKNOWN_CHECKSUM)
    assert sf.getChecksum() == b"abcde000"

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
                pyopenms.TransformationModelBSpline,
                pyopenms.TransformationModelInterpolated]:
        p = pyopenms.Param()
        data = [ pyopenms.TM_DataPoint(9.0, 8.9),
                 pyopenms.TM_DataPoint(5.0, 6.0),
                 pyopenms.TM_DataPoint(8.0, 8.0) ]
        mod = clz(data, p)
        mod.evaluate(7.0)
        mod.getDefaultParameters(p)

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
    fh.store(b"test.transformationXML", td)
    fh.load(b"test.transformationXML", td, True)
    assert td.getDataPoints() == []

@report
def testIBSpectraFile():
    """
    @tests:
     IBSpectraFile.__init__
     IBSpectraFile.store
    """
    fh = pyopenms.IBSpectraFile()
    cmap = pyopenms.ConsensusMap()
    correctError = False
    try:
        fh.store( pyopenms.String(b"test.ibspectra.file"), cmap)
        assert False
    except RuntimeError:
        correctError = True

    assert correctError 

@report
def testSwathFile():
    """
    @tests:
     SwathFile.__init__
     SwathFile.store
    """
    fh = pyopenms.SwathFile()

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
    assert isinstance( pyopenms.VersionInfo.getVersion(), bytes)
    assert isinstance( pyopenms.VersionInfo.getRevision(), bytes)
    assert isinstance( pyopenms.VersionInfo.getTime(), bytes)

    vd = pyopenms.VersionDetails.create(b"19.2.1")
    assert vd.version_major == 19
    assert vd.version_minor == 2
    assert vd.version_patch == 1

    assert vd == vd
    assert not vd < vd
    assert not vd > vd

    assert  isinstance(pyopenms.version.version, str)

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

    inst.name = b"test"
    inst.value = b"test"
    inst.cvRef = b"test"
    inst.cvAcc = b"test"
    inst.unitRef = b"test"
    inst.unitAcc = b"test"
    inst.binary = b"test"
    inst.qualityRef = b"test"
    inst.colTypes = [ b"test", b"test2"]
    inst.tableRows = [ [b"test", b"test2"], [b"otherTest"] ]

    assert inst.tableRows[1][0] == b"otherTest"

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
def testKernelMassTrace():
    trace = pyopenms.Kernel_MassTrace()

    assert trace.getSize is not None
    assert trace.getLabel is not None
    assert trace.setLabel is not None

    assert trace.getCentroidMZ is not None
    assert trace.getCentroidRT is not None
    assert trace.getCentroidSD is not None
    assert trace.getFWHM is not None
    assert trace.getTraceLength is not None
    assert trace.getFWHMborders is not None
    assert trace.getSmoothedIntensities is not None
    assert trace.getAverageMS1CycleTime is not None

    assert trace.computeSmoothedPeakArea is not None
    assert trace.computePeakArea is not None
    assert trace.findMaxByIntPeak is not None
    assert trace.estimateFWHM is not None
    assert trace.computeFwhmArea is not None
    assert trace.computeFwhmAreaSmooth is not None
    # assert trace.computeFwhmAreaRobust is not None
    # assert trace.computeFwhmAreaSmoothRobust is not None
    assert trace.getIntensity is not None
    assert trace.getMaxIntensity is not None

    assert trace.getConvexhull is not None

    assert trace.setCentroidSD is not None
    assert trace.setSmoothedIntensities is not None
    assert trace.updateSmoothedMaxRT is not None
    assert trace.updateWeightedMeanRT is not None
    assert trace.updateSmoothedWeightedMeanRT is not None
    assert trace.updateMedianRT is not None
    assert trace.updateMedianMZ is not None
    assert trace.updateMeanMZ is not None
    assert trace.updateWeightedMeanMZ is not None
    assert trace.updateWeightedMZsd is not None

    s = trace.getSize()

@report
def testElutionPeakDetection():
    detection = pyopenms.ElutionPeakDetection()

    assert detection.detectPeaks is not None
    assert detection.filterByPeakWidth is not None
    assert detection.computeMassTraceNoise is not None
    assert detection.computeMassTraceSNR is not None
    assert detection.computeApexSNR is not None
    assert detection.findLocalExtrema  is not None
    assert detection.smoothData  is not None

    trace = pyopenms.Kernel_MassTrace()
    detection.smoothData(trace, 4)

@report
def testIndexedMzMLDecoder():
    decoder = pyopenms.IndexedMzMLDecoder()

    try:
        pos = decoder.findIndexListOffset(b"abcde", 100)
        raise Exception("Should raise an error")
    except RuntimeError:
        pass

def test_streampos():
    p = long(pyopenms.streampos())
    assert isinstance(p, long), "got %r" % p

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

def test_BSpline2d():

    x = [1.0, 6.0, 8.0, 10.0, 15.0]
    y = [2.0, 5.0, 6.0, 12.0, 13.0]
    spline = pyopenms.BSpline2d(x,y,0, pyopenms.BoundaryCondition.BC_ZERO_ENDPOINTS, 0)

    assert spline.ok()
    assert abs(spline.eval(6.0) - 5.0 < 0.01)
    assert abs(spline.derivative(6.0) - 5.0 < 0.01)

    y_new = [4.0, 5.0, 6.0, 12.0, 13.0]
    spline.solve(y_new)

    assert spline.ok()
    assert abs(spline.eval(6.0) - 5.0 < 0.01)


@report
def testConsensusIDAlgorithmAverage():
    algo = pyopenms.ConsensusIDAlgorithmAverage()
    assert algo.apply

@report
def testConsensusIDAlgorithmBest():
    algo = pyopenms.ConsensusIDAlgorithmBest()
    assert algo.apply

@report
def testConsensusIDAlgorithmIdentity():
    algo = pyopenms.ConsensusIDAlgorithmIdentity()
    assert algo.apply

@report
def testConsensusIDAlgorithmPEPIons():
    algo = pyopenms.ConsensusIDAlgorithmPEPIons()
    assert algo.apply

@report
def testConsensusIDAlgorithmPEPMatrix():
    algo = pyopenms.ConsensusIDAlgorithmPEPMatrix()
    assert algo.apply

@report
def testConsensusIDAlgorithmRanks():
    algo = pyopenms.ConsensusIDAlgorithmRanks()
    assert algo.apply

@report
def testConsensusIDAlgorithmSimilarity():
    algo = pyopenms.ConsensusIDAlgorithmSimilarity()
    assert algo.apply

@report
def testConsensusIDAlgorithmWorst():
    algo = pyopenms.ConsensusIDAlgorithmWorst()
    assert algo.apply

@report
def testDigestionEnzymeProtein():
    f = pyopenms.EmpiricalFormula()

    regex_description = b""
    psi_id = b""
    xtandem_id = b""
    comet_id = 0
    omssa_id = 0
    e = pyopenms.DigestionEnzymeProtein(b"testEnzyme", "K", set([]), regex_description,
                                 f, f, psi_id, xtandem_id, comet_id, omssa_id)

@report
def testMRMAssay():
    e = pyopenms.MRMAssay()
    assert e

@report
def testMRMIonSeries():
    e = pyopenms.MRMIonSeries()
    assert e

@report
def testPeptideIndexing():
    e = pyopenms.PeptideIndexing()
    assert e

@report
def testPeptideProteinResolution():
    e = pyopenms.PeptideProteinResolution(False)
    assert e

@report
def testPercolatorOutfile():
    e = pyopenms.PercolatorOutfile()
    assert e



@report
def testHiddenMarkovModel():
    hmm = pyopenms.HiddenMarkovModel()
    assert hmm

    assert hmm.getNumberOfStates() == 0

    ss = s("testState")
    hmm.addNewState(ss)

    assert hmm.getNumberOfStates() == 1

    e = pyopenms.HMMState()
    # hmm.addNewState(e) # Segfault !

    r = hmm.getState(s("testState"))
    assert r
    ## assert r == ss # requires ==

@report
def testHMMState():
    e = pyopenms.HMMState()
    assert e
    e.setName(s("somename"))
    assert e.getName() == b"somename", e.getName()
    e.setHidden(True)
    assert e.isHidden()

    pre = pyopenms.HMMState()
    pre.setName(s("pre"))
    suc = pyopenms.HMMState()
    suc.setName(s("suc"))

    e.addPredecessorState(pre)
    e.addSuccessorState(suc)

    assert e.getPredecessorStates()
    assert e.getSuccessorStates()



@report
def testProteaseDB():
    edb = pyopenms.ProteaseDB()

    f = pyopenms.EmpiricalFormula()
    synonyms = set([b"dummy", b"other"])

    assert edb.hasEnzyme(pyopenms.String("Trypsin"))

    trypsin = edb.getEnzyme(pyopenms.String("Trypsin"))

    names = []
    edb.getAllNames(names)
    assert b"Trypsin" in names


@report
def testElementDB():
    edb = pyopenms.ElementDB()
    del edb

    # create a second instance of ElementDB without anything bad happening
    edb = pyopenms.ElementDB()

    assert edb.hasElement(16)
    edb.hasElement(pyopenms.String(b"O"))

    e = edb.getElement(16)

    assert e.getName() == b"Sulfur"
    assert e.getSymbol() == b"S"
    assert e.getIsotopeDistribution()

    e2 = edb.getElement(pyopenms.String(b"O"))

    assert e2.getName() == b"Oxygen"
    assert e2.getSymbol() == b"O"
    assert e2.getIsotopeDistribution()

    # assert e == e2

    #  not yet implemented
    #
    # const Map[ String, Element * ]  getNames() nogil except +
    # const Map[ String, Element * ] getSymbols() nogil except +
    # const Map[unsigned int, Element * ] getAtomicNumbers() nogil except +


@report
def testDPosition():
    dp = pyopenms.DPosition1()
    dp = pyopenms.DPosition1(1.0)
    assert dp[0] == 1.0

    dp = pyopenms.DPosition2()
    dp = pyopenms.DPosition2(1.0, 2.0)

    assert dp[0] == 1.0
    assert dp[1] == 2.0

@report
def testResidueDB():
    rdb = pyopenms.ResidueDB()
    del rdb

    # create a second instance of ResidueDB without anything bad happening
    rdb = pyopenms.ResidueDB()

    assert rdb.getNumberOfResidues() >= 20
    assert len(rdb.getResidueSets() ) >= 1
    el = rdb.getResidues(pyopenms.String(rdb.getResidueSets().pop()))

    assert len(el) >= 1

    assert rdb.hasResidue(s(b"Glycine"))
    glycine = rdb.getResidue(s(b"Glycine"))

    nrr = rdb.getNumberOfResidues()

    r = pyopenms.Residue()
    rdb.addResidue(r)
    assert rdb.getNumberOfResidues() == nrr+1

@report
def testModificationsDB():
    mdb = pyopenms.ModificationsDB()
    del mdb

    # create a second instance of ModificationsDB without anything bad happening
    mdb = pyopenms.ModificationsDB()

    assert mdb.getNumberOfModifications() > 1
    m = mdb.getModification(1)

    assert mdb.getNumberOfModifications() > 1
    m = mdb.getModification(1)
    assert m is not None

    mods = set([])
    mdb.searchModifications(mods, s("Phosphorylation"), s("T"), pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert len(mods) == 1

    mods = set([])
    mdb.searchModifications(mods, s("NIC"), s("T"), pyopenms.ResidueModification.TermSpecificity.N_TERM)
    assert len(mods) == 1

    mods = set([])
    mdb.searchModifications(mods, s("NIC"), s("T"), pyopenms.ResidueModification.TermSpecificity.N_TERM)
    assert len(mods) == 1

    mods = set([])
    mdb.searchModifications(mods, s("Acetyl"), s("T"), pyopenms.ResidueModification.TermSpecificity.N_TERM)
    assert len(mods) == 1
    assert list(mods)[0].getFullId() == b"Acetyl (N-term)"

    m = mdb.getModification(s("Carboxymethyl (C)"), "", pyopenms.ResidueModification.TermSpecificity.NUMBER_OF_TERM_SPECIFICITY)
    assert m.getFullId() == b"Carboxymethyl (C)"

    m = mdb.getModification( s("Phosphorylation"), s("S"), pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m.getId() == b"Phospho"

    # get out all mods (there should be many, some known ones as well!)
    mods = []
    m = mdb.getAllSearchModifications(mods)
    assert len(mods) > 100

    assert b"Phospho (S)" in mods
    assert b"Sulfo (S)" in mods
    assert not (b"Phospho" in mods)

    # search for specific modifications by mass
    m = mdb.getBestModificationByDiffMonoMass( 80.0, 1.0, b"T", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getId() == b"Phospho"
    assert m.getFullName() == b"Phosphorylation"
    assert m.getUniModAccession() == b"UniMod:21"

    m = mdb.getBestModificationByDiffMonoMass(80, 100, b"T", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getId() == b"Phospho"
    assert m.getFullName() == b"Phosphorylation"
    assert m.getUniModAccession() == b"UniMod:21"

    m = mdb.getBestModificationByDiffMonoMass(16, 1.0, b"M", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getId() == b"Oxidation", m.getId()
    assert m.getFullName() == b"Oxidation or Hydroxylation", m.getFullName()
    assert m.getUniModAccession() == b"UniMod:35"

    ###

    m = mdb.getBestModificationByMonoMass(80, 20, b"T", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getId() == b"MOD:00439"
    assert m.getFullName() == b"O-phospho-L-threonine with neutral loss of phosphate", m.getFullName() # something crazy
    assert m.getUniModAccession() == b"" # no unimod for crazyness ...

    m = mdb.getBestModificationByMonoMass(147, 20, b"M", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getUniModAccession() == b"", m.getUniModAccession()
    assert m.getId() == b"MOD:00719", m.getId()
    assert m.getFullName() == b"oxidation to L-methionine sulfoxide", m.getFullName()

    m = mdb.getBestModificationByMonoMass( 96, 20, b"T", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getId() == b"MOD:00252", m.getId()
    assert m.getFullName() == b"keratan sulfate D-glucuronosyl-D-galactosyl-D-galactosyl-D-xylosyl-L-threonine", m.getFullName() # something crazy
    assert m.getUniModAccession() == b"", m.getUniModAccession() # no unimod for crazyness ...

    # Test NULL ptr
    m = mdb.getBestModificationByMonoMass( 999999999, 0.20, b"T", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is None

def testString():
    pystr = pyopenms.String()
    pystr = pyopenms.String("blah")
    assert (pystr.toString() == "blah")
    pystr = pyopenms.String(b"blah")
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
    # pystr = pyopenms.String(b"bläh") # Not Python 3 compatible
    pystr = pyopenms.String(u"bläh")
    assert (pystr.toString() == u"bläh")
    pystr = pyopenms.String(pystr)
    assert (pystr.toString() == u"bläh")
    cstr = pystr.c_str()

    # Printing should work ...
    print(cstr)
    print(pystr)
    print(pystr.toString().encode("utf8"))

    assert len(pystr.toString()) == 4
    assert len(pystr.c_str()) == 5 # C does not know about Unicode, so be careful with c_str
    print(pystr) # this prints the C string, due to Py 2/3 compatibility
    print(pystr.toString().encode("utf8")) # this prints the correct String

    pystr1 = pyopenms.String("bläh")
    pystr2 = pyopenms.String("bläh")
    assert(pystr1 == pystr2)

    pystr1 = pyopenms.String(u"bläh")
    pystr2 = pyopenms.String(u"bläh")
    assert(pystr1 == pystr2)


