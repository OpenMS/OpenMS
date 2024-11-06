#!/usr/bin/env python
# -*- coding: utf-8  -*-
from __future__ import print_function

import pyopenms
import copy
import os
import logging

from pyopenms import String as s
import numpy as np
import pandas as pd

print("IMPORTED ", pyopenms.__file__)

try:
    long
except NameError:
    long = int

from functools import wraps

import sys
def _testStrOutput(input_str):
    if sys.version_info[0] < 3:
        assert isinstance(input_str, unicode)
    else:
        assert isinstance( input_str, str)

def report(f):
    @wraps(f)
    def wrapper(*a, **kw):
        print("run ", f.__name__)
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
    assert len(keys) and all(isinstance(k, bytes) for k in keys)
    assert what.getMetaValue(keys[0]) == 42

    assert what.metaValueExists("key")
    what.removeMetaValue("key")

    keys = []
    what.getKeys(keys)
    assert what.getMetaValue(keys[0]) == 42

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
    @tests: ProgressLogger
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
    @tests: SpectrumAlignment
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
    @tests: AASequence
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
    aas = pyopenms.AASequence.fromString("DFPIANGER")
    assert aas.getCTerminalModificationName() == ""
    assert aas.getNTerminalModificationName() == ""
    aas.setCTerminalModification("")
    aas.setNTerminalModification("")
    assert aas.toString() == "DFPIANGER"
    assert aas.toUnmodifiedString() == "DFPIANGER"
    aas = pyopenms.AASequence.fromStringPermissive("DFPIANGER", True)
    assert aas.toString() == "DFPIANGER"
    assert aas.toUnmodifiedString() == "DFPIANGER"

    seq = pyopenms.AASequence.fromString("PEPTIDESEKUEM(Oxidation)CER")
    assert seq.toString() == "PEPTIDESEKUEM(Oxidation)CER"
    assert seq.toUnmodifiedString() == "PEPTIDESEKUEMCER"
    assert seq.toBracketString() == "PEPTIDESEKUEM[147]CER"
    assert seq.toBracketString(True) == "PEPTIDESEKUEM[147]CER"

    assert seq.toBracketString(False) == "PEPTIDESEKUEM[147.03540001709996]CER" or \
           seq.toBracketString(False) == "PEPTIDESEKUEM[147.035400017100017]CER"

    assert seq.toBracketString(False) == "PEPTIDESEKUEM[147.03540001709996]CER" or \
           seq.toBracketString(False) == "PEPTIDESEKUEM[147.035400017100017]CER"

    assert seq.toUniModString() == "PEPTIDESEKUEM(UniMod:35)CER"
    assert seq.isModified()
    assert not seq.hasCTerminalModification()
    assert not seq.hasNTerminalModification()
    assert not seq.empty()

    # has selenocysteine
    assert seq.getResidue(1) is not None
    assert seq.size() == 16

    # test exception forwarding from C++ to python
    # classes derived from std::runtime_exception can be caught in python
    try:
        seq.getResidue(1000) # does not exist
    except RuntimeError:
        print("Exception successfully triggered.")
    else:
        print("Error: Exception not triggered.")
        assert False
    assert seq.getFormula(pyopenms.Residue.ResidueType.Full, 0) == pyopenms.EmpiricalFormula("C75H122N20O32S2Se1")
    assert abs(seq.getMonoWeight(pyopenms.Residue.ResidueType.Full, 0) - 1958.7140766518) < 1e-5
    # assert seq.has(pyopenms.ResidueDB.getResidue("P"))

    
@report
def testElement():
    """
    @tests: Element
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

    e = pyopenms.Element()
    e.setSymbol("blah")
    e.setSymbol("blah")
    e.setSymbol(u"blah")
    e.setSymbol(str("blah"))
    oms_string = s("blu")
    e.setSymbol(oms_string)
    assert oms_string
    assert oms_string.toString() == "blu"

    evil = u"blü"
    evil8 = evil.encode("utf8")
    evil1 = evil.encode("latin1")


    e.setSymbol(evil.encode("utf8"))
    assert e.getSymbol() == u"blü"
    e.setSymbol(evil.encode("latin1"))
    assert e.getSymbol().decode("latin1") == u"blü"

    # If we get the raw symbols, we get bytes (which we would need to decode first)
    e.setSymbol(evil8.decode("utf8"))
    # assert e.getSymbol() == 'bl\xc3\xbc', e.getSymbol()
    assert e.getSymbol() == u"blü" #.encode("utf8")
    # OpenMS strings, however, understand the decoding
    assert s(e.getSymbol()) == s(u"blü")
    assert s(e.getSymbol()).toString() == u"blü"

    # What if you use the wrong decoding ?
    e.setSymbol(evil1)
    assert e.getSymbol().decode("latin1") == u"blü"
    e.setSymbol(evil8)
    assert e.getSymbol() == u"blü"

@report
def testResidue():
    """
    @tests: Residue
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
    @tests: IsotopeDistribution
     IsotopeDistribution.__init__
    """
    ins = pyopenms.IsotopeDistribution()

    ins.getMax()
    ins.getMin()
    ins.size()
    ins.clear()
    ins.renormalize()
    ins.trimLeft(6.0)
    ins.trimRight(8.0)

    ins.clear()
    ins.insert(1, 2)
    ins.insert(6, 5)

    assert ins.size() == 2

    for p in ins:
        print(p)

@report
def testFineIsotopePatternGenerator():
    """
    @tests: FineIsotopePatternGenerator
    """

    iso = pyopenms.FineIsotopePatternGenerator()
    iso.setThreshold(1e-5)
    iso.setAbsolute(True)
    assert iso.getAbsolute() 

    methanol = pyopenms.EmpiricalFormula("CH3OH")
    water = pyopenms.EmpiricalFormula("H2O")
    mw = methanol + water
    iso_dist = mw.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-20, False, False))
    assert len(iso_dist.getContainer()) == 56
    iso_dist = mw.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-200, False, False))
    assert len(iso_dist.getContainer()) == 84

    c100 = pyopenms.EmpiricalFormula("C100")
    iso_dist = c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-200, False, False))
    assert len(iso_dist.getContainer()) == 101
    assert c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-2, False, False)).size() == 6
    assert c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-2, False, True)).size() == 5
    assert c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-2, True, False)).size() == 5
    assert c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-2, True, True)).size() == 5

    assert c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-10, False, False)).size() == 14
    assert c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-10, False, True)).size() == 13
    assert c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-10, True, False)).size() == 10
    assert c100.getIsotopeDistribution(pyopenms.FineIsotopePatternGenerator(1e-10, True, True)).size() == 10

    iso = pyopenms.FineIsotopePatternGenerator(1e-5, False, False)
    isod = iso.run(methanol)
    assert len(isod.getContainer()) == 6
    assert abs(isod.getContainer()[0].getMZ() - 32.0262151276) < 1e-5
    assert isod.getContainer()[0].getIntensity() - 0.986442089081 < 1e-5

@report
def testCoarseIsotopePatternGenerator():
    """
    @tests: CoarseIsotopePatternGenerator
    CoarseIsotopePatternGenerator.__init__
    CoarseIsotopePatternGenerator.getMaxIsotope()
    CoarseIsotopePatternGenerator.setMaxIsotope()
    CoarseIsotopePatternGenerator.estimateFromPeptideWeight()
    """

    iso = pyopenms.CoarseIsotopePatternGenerator()
    iso.setMaxIsotope(5)
    assert iso.getMaxIsotope() == 5
    res = iso.estimateFromPeptideWeight(500)

    methanol = pyopenms.EmpiricalFormula("CH3OH")
    water = pyopenms.EmpiricalFormula("H2O")
    mw = methanol + water
    iso_dist = mw.getIsotopeDistribution(pyopenms.CoarseIsotopePatternGenerator(3))
    assert len(iso_dist.getContainer()) == 3, len(iso_dist.getContainer())
    iso_dist = mw.getIsotopeDistribution(pyopenms.CoarseIsotopePatternGenerator(0))
    assert len(iso_dist.getContainer()) == 18, len(iso_dist.getContainer()) 

    iso = pyopenms.CoarseIsotopePatternGenerator(10)
    isod = iso.run(methanol)
    assert len(isod.getContainer()) == 10, len(isod.getContainer()) 
    assert abs(isod.getContainer()[0].getMZ() - 32.0262151276) < 1e-5
    assert isod.getContainer()[0].getIntensity() - 0.986442089081 < 1e-5
    
@report
def testEmpiricalFormula():
    """
    @tests: EmpiricalFormula 
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
    ins.getIsotopeDistribution(pyopenms.CoarseIsotopePatternGenerator(0))
    # ins.getNumberOf(0)
    # ins.getNumberOf("test")
    ins.getNumberOfAtoms()
    ins.setCharge(2)
    ins.getCharge()
    ins.toString()
    ins.isEmpty()
    ins.isCharged()
    ins.hasElement( pyopenms.Element() )

    ef = pyopenms.EmpiricalFormula("C2H5")
    s = ef.toString()
    assert s == "C2H5"
    m = ef.getElementalComposition()
    assert m[b"C"] == 2
    assert m[b"H"] == 5
    assert ef.getNumberOfAtoms() == 7

@report
def testModificationDefinitionsSet():
    """
    @tests: ModificationDefinitionsSet
     ModificationDefinitionsSet.__init__
    """
    empty = pyopenms.ModificationDefinitionsSet()
    fixed = [b"Carbamidomethyl"]
    variable = [b"Oxidation"]
    full = pyopenms.ModificationDefinitionsSet(fixed, variable)

@report
def test_AcquisitionInfo():
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
def test_BaseFeature():
    """
    @tests: BaseFeature
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
    @tests: AnnotationState
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
    @tests: ChecksumType
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
    @tests: ChromatogramPeak
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
    @tests: ChromatogramTools
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
    @tests: ConsensusFeature
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
    @tests: ConsensusMap
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
    
    f = pyopenms.ConsensusFeature()
    f.setCharge(1)
    f.setQuality(2.0)
    f.setWidth(4.0)
    m.push_back(f)
    m.push_back(f)
    intydf = m.get_intensity_df()
    metadf = m.get_metadata_df()
    assert intydf.shape[0] == 2
    assert metadf.shape[0] == 2

    assert m == m
    assert not m != m

@report
def testConsensusXMLFile():
    """
    @tests: ConsensusXMLFile
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
    @tests: XTandemXMLFile
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


    f.setTaxon("testTaxon")
    assert f.getTaxon() == "testTaxon"

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
    @tests: SignalToNoiseEstimatorMedian
     SignalToNoiseEstimatorMedian.__init__
    """
    f = pyopenms.SignalToNoiseEstimatorMedian()
    assert f.init is not None
    assert f.getSignalToNoise is not None

@report
def testSignalToNoiseEstimatorMedianChrom():
    """
    @tests: SignalToNoiseEstimatorMedianChrom
     SignalToNoiseEstimatorMedianChrom.__init__
    """
    f = pyopenms.SignalToNoiseEstimatorMedianChrom()
    assert f.init is not None
    assert f.getSignalToNoise is not None

@report
def testConvexHull2D():
    """
    @tests: ConvexHull2D
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

    _testMetaInfoInterface(dp)

    assert dp == dp
    assert not dp != dp

    # assert isinstance(dp.getCompletionTime().getDate(), bytes)
    # assert isinstance(dp.getCompletionTime().getTime(), bytes)
    dp.clearMetaInfo()
    k = []
    dp.getKeys(k)
    assert k == []
    dp.getMetaValue
    ac = dp.getProcessingActions()
    assert ac == set(())
    ac = set([ pyopenms.DataProcessing.ProcessingAction.PEAK_PICKING, pyopenms.DataProcessing.ProcessingAction.BASELINE_REDUCTION])
    dp.setProcessingActions(ac)
    assert len(dp.getProcessingActions() ) == 2
    _testStrOutput(dp.getSoftware().getName())
    _testStrOutput(dp.getSoftware().getVersion())
    dp.isMetaEmpty()
    dp.metaValueExists
    dp.removeMetaValue
    # dp.setCompletionTime(pyopenms.DateTime.now())
    s = dp.getSoftware()
    s.setName("pyopenms")
    dp.setSoftware(s)

    assert dp.getSoftware().getName() == "pyopenms"


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
def testAdduct():
    """
    @tests: Adduct
     Adduct.__init__
    """
    a = pyopenms.Adduct()

@report
def testGaussFitter():
    """
    @tests: GaussFitter
     GaussFitter.__init__
    """
    ins = pyopenms.GaussFitter()

@report
def testGaussFitResult():
    """
    @tests: GaussFitResult
     GaussFitResult.__init__
    """
    ins = pyopenms.GaussFitResult(0.0, 0.0, 0.0)
    ins.A = 5.0
    ins.x0 = 5.0
    ins.sigma = 5.0

@report
def testChargePair():
    """
    @tests: ChargePair
     ChargePair.__init__
    """
    a = pyopenms.ChargePair()

@report
def testCompomer():
    """
    @tests: Compomer
     Compomer.__init__
    """
    a = pyopenms.Compomer()

@report
def testCVMappings():
    """
    @tests: CVMappings
     CVMappings.__init__
    """
    val = pyopenms.CVMappings()

@report
def testCVMappingFile():
    """
    @tests: CVMappingFile
     CVMappingFile.__init__
    """
    val = pyopenms.CVMappingFile()

    assert pyopenms.CVMappingFile().load

@report
def testControlledVocabulary():
    """
    @tests: ControlledVocabulary
     ControlledVocabulary.__init__
    """
    val = pyopenms.ControlledVocabulary()

    assert pyopenms.ControlledVocabulary().loadFromOBO

@report
def testSemanticValidator():
    """
    @tests: SemanticValidator
     SemanticValidator.__init__
    """
    m = pyopenms.CVMappings()
    cv = pyopenms.ControlledVocabulary()

    val = pyopenms.SemanticValidator(m, cv)

    assert val.validate is not None
    assert val.setCheckTermValueTypes is not None
    assert val.setCheckUnits is not None


# @report
# def testDateTime():
#     """
#     @tests: DateTime
#      DateTime.__init__
#      DateTime.getDate
#      DateTime.getTime
#      DateTime.now
#     """
#     d = pyopenms.DateTime()
#     assert isinstance( d.getDate(), bytes)
#     assert isinstance( d.getTime(), bytes)
#     d = pyopenms.DateTime.now()
#     assert isinstance( d.getDate(), bytes)
#     assert isinstance( d.getTime(), bytes)
# 
#     d.clear()
#     d.set("01.01.2001 11:11:11")
#     assert d.get() == "2001-01-01 11:11:11"

@report
def testFeature():
    """
    @tests: Feature
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
def testFeatureFileOptions():
    """
    @tests: FeatureFileOptions
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
    @tests: Param
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
     Param.getValidStrings
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
        # only set the section description if there are actually two or more sections
        if len(k.split(b":")) < 2: continue
        f = k.split(b":")[0]
        p.setSectionDescription(f, k)
        # TODO: keys inside maps are not yet properly decoded
        assert p.getSectionDescription(f) == k.decode()

        assert p.get(k) is not None

    assert len(p.values()) == len([p[k] for k in p.keys()])
    assert sorted(p.items()) == sorted((k, p[k]) for k in p.keys())

    assert not p.exists("asdflkj01231321321v")
    p.addTag(k, "a")
    p.addTags(k, [b"", b"c"])
    assert sorted(p.getTags(k)) == [b"", b"a", b"c"]
    p.clearTags(k)
    assert p.getTags(k) == []

    pn = pyopenms.Param()
    pn.insert("master:", p)
    assert pn.exists(b"master:"+k)

    p1 = pn.copy("master:", True)
    assert p1 == p

    p1.update(p)
    p1.update(p,0)
    p1.update(p,1)
    p1.update(dd)

    p.setValidStrings
    p.getValidStrings
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

    assert p1.get(b"abcde", 7) == 7




@report
def testFeatureFinderAlgorithmPicked():
    """
    @tests: FeatureFinderAlgorithmPicked
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

    assert ff.getName() == "FeatureFinderAlgorithmPicked"

    ff.setParameters(pyopenms.Param())

    ff.setName("test")
    assert ff.getName() == "test"


@report
def testExperimentalSettings():
    """
    @tests: ExperimentalSettings
     ExperimentalSettings.__init__
    """
    ff = pyopenms.ExperimentalSettings()

@report
def testFeatureDeconvolution():
    """
    @tests: FeatureDeconvolution
     FeatureDeconvolution.__init__
    """
    ff = pyopenms.FeatureDeconvolution()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.FeatureDeconvolution().compute is not None

@report
def testInternalCalibration():
    """
    @tests: InternalCalibration
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
    @tests: testItraqConstants
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
    @tests: LinearResampler
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
    @tests: PeptideAndProteinQuant
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
    @tests: SeedListGenerator
     SeedListGenerator.__init__
    """
    ff = pyopenms.SeedListGenerator()
    p = ff.getDefaults()
    _testParam(p)

    # TODO 
    # assert pyopenms.SeedListGenerator().compute is not None

# TODO: re-enable as soon as ConsensusIDAlgorithm classes are wrapped
# @report
# def testConsensusID():
#     """
#     @tests: ConsensusID
#      ConsensusID.__init__
#     """
#     ff = pyopenms.ConsensusID()
#     p = ff.getDefaults()
#     _testParam(p)

#     assert pyopenms.ConsensusID().apply is not None

@report
def testFalseDiscoveryRate():
    """
    @tests: FalseDiscoveryRate
     FalseDiscoveryRate.__init__
    """
    ff = pyopenms.FalseDiscoveryRate()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.FalseDiscoveryRate().apply is not None

@report
def testIDFilter():
    """
    @tests: IDFilter
     IDFilter.__init__
    """
    ff = pyopenms.IDFilter()

    # assert pyopenms.IDFilter().apply is not None

@report
def testPosteriorErrorProbabilityModel():
    """
    @tests: PosteriorErrorProbabilityModel
     PosteriorErrorProbabilityModel.__init__
    """
    model = pyopenms.PosteriorErrorProbabilityModel()
    p = model.getDefaults()
    _testParam(p)

    assert pyopenms.PosteriorErrorProbabilityModel().fit is not None
    assert pyopenms.PosteriorErrorProbabilityModel().computeProbability is not None

    scores = [float(i) for i in range(10)]
    model.fit(scores, "none")
    model.fit(scores, scores, "none")

    model.fillLogDensities(scores, scores, scores)

    assert model.computeLogLikelihood is not None
    assert model.pos_neg_mean_weighted_posteriors is not None

    GaussFitResult = model.getCorrectlyAssignedFitResult()
    GaussFitResult = model.getIncorrectlyAssignedFitResult()
    model.getNegativePrior()
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
    @tests: SeedListGenerator
     SeedListGenerator.__init__
    """
    ff = pyopenms.SeedListGenerator()

    # TODO 
    # assert pyopenms.SeedListGenerator().generateSeedList is not None

@report
def testConsensusMapNormalizerAlgorithmMedian():
    """
    @tests: ConsensusMapNormalizerAlgorithmMedian
     ConsensusMapNormalizerAlgorithmMedian.__init__
    """
    ff = pyopenms.ConsensusMapNormalizerAlgorithmMedian()

    assert pyopenms.ConsensusMapNormalizerAlgorithmMedian().normalizeMaps is not None

@report
def testConsensusMapNormalizerAlgorithmQuantile():
    """
    @tests: ConsensusMapNormalizerAlgorithmQuantile
     ConsensusMapNormalizerAlgorithmQuantile.__init__
    """
    ff = pyopenms.ConsensusMapNormalizerAlgorithmQuantile()

    assert pyopenms.ConsensusMapNormalizerAlgorithmQuantile().normalizeMaps is not None

@report
def testConsensusMapNormalizerAlgorithmThreshold():
    """
    @tests: ConsensusMapNormalizerAlgorithmThreshold
     ConsensusMapNormalizerAlgorithmThreshold.__init__
    """
    ff = pyopenms.ConsensusMapNormalizerAlgorithmThreshold()

    assert pyopenms.ConsensusMapNormalizerAlgorithmThreshold().computeCorrelation is not None
    assert pyopenms.ConsensusMapNormalizerAlgorithmThreshold().normalizeMaps is not None


@report
def testAScore():
    """
    @tests: AScore
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
    @tests: IDRipper
     IDRipper.__init__
     IDRipper.rip
    """
    ff = pyopenms.IDRipper()

    assert pyopenms.IDRipper().rip is not None

@report
def testFASTAFile():
    """
    @tests: FASTAFile
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
    @tests: FASTAEntry
     FASTAEntry.__init__
    """
    ff = pyopenms.FASTAEntry()

@report
def testInternalCalibration():
    """
    @tests: InternalCalibration
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
    @tests: ProteaseDigestion
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
def testIDDecoyProbability():
    """
    @tests: IDDecoyProbability
      IDDecoyProbability.__init__
    """
    ff = pyopenms.IDDecoyProbability()

    assert pyopenms.IDDecoyProbability().apply is not None

@report
def testFeatureGrouping():
    """
    @tests: FeatureGroupingAlgorithm
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
    @tests: FeatureMap
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

    assert isinstance(fm2.getMinRT(), float)
    assert isinstance(fm2.getMinRT(), float)
    assert isinstance(fm2.getMaxMZ(), float)
    assert isinstance(fm2.getMaxMZ(), float)
    assert isinstance(fm2.getMinIntensity(), float)
    assert isinstance(fm2.getMaxIntensity(), float)

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
    @tests: FeatureXMLFile
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

    f = pyopenms.Feature()
    f.setMZ(200)
    f.setCharge(1)
    f.setRT(10)
    f.setIntensity(10000)
    f.setOverallQuality(10)

    ch = pyopenms.ConvexHull2D()
    ch.setHullPoints(np.asarray([[8,199],[12,201]], dtype='f'))
    f.setConvexHulls([ch])

    f.setMetaValue(b'mv1', 1)
    f.setMetaValue(b'mv2', 2)

    f.setMetaValue('spectrum_native_id', 'spectrum=123')
    pep_id = pyopenms.PeptideIdentification()
    pep_id.insertHit(pyopenms.PeptideHit())
    f.setPeptideIdentifications([pep_id])

    fm.push_back(f)

    f.setMetaValue('spectrum_native_id', 'spectrum=124')
    fm.push_back(f)

    assert len(fm.get_assigned_peptide_identifications()) == 2
    assert fm.get_df(meta_values='all').shape == (2, 16)
    assert fm.get_df(meta_values='all', export_peptide_identifications=False).shape == (2, 12)

    assert pd.merge(fm.get_df(), pyopenms.peptide_identifications_to_df(fm.get_assigned_peptide_identifications()),
                on = ['feature_id', 'ID_native_id', 'ID_filename']).shape == (2,24)

    fm = pyopenms.FeatureMap()
    pyopenms.FeatureXMLFile().load(os.path.join(os.environ['OPENMS_DATA_PATH'], 'examples/FRACTIONS/BSA1_F1_idmapped.featureXML'), fm)

    assert pd.merge(fm.get_df(), pyopenms.peptide_identifications_to_df(fm.get_assigned_peptide_identifications()),
                    on = ['feature_id', 'ID_native_id', 'ID_filename']).shape == (15,26)

    fh = pyopenms.FeatureXMLFile()
    fh.store("test.featureXML", fm)
    fh.load("test.featureXML", fm)

    fh = pyopenms.FileHandler()
    fh.loadFeatures("test.featureXML", fm)

@report
def testFileDescription():
    """
    @tests: ColumnHeader
     ColumnHeader.__init__
     ColumnHeader.filename
     ColumnHeader.label
     ColumnHeader.size
     ColumnHeader.unique_id
    """
    fd = pyopenms.ColumnHeader()
    _testStrOutput(fd.filename)
    _testStrOutput(fd.label)
    assert isinstance(fd.size, int)
    # assert isinstance(fd.unique_id, (long, int, bytes))

@report
def testFileHandler():
    """
    @tests: FileHandler
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
def testCachedMzML():
    """
    """
    mse = pyopenms.MSExperiment()
    s = pyopenms.MSSpectrum()
    mse.addSpectrum(s)

    # First load data and cache to disk
    pyopenms.CachedmzML.store("myCache.mzML", mse)

    # Now load data
    cfile = pyopenms.CachedmzML()
    pyopenms.CachedmzML.load("myCache.mzML", cfile)

    meta_data = cfile.getMetaData()
    assert cfile.getNrChromatograms() ==0
    assert cfile.getNrSpectra() == 1

@report
def testIndexedMzMLFile():
    """
    """
    mse = pyopenms.MSExperiment()
    s = pyopenms.MSSpectrum()
    mse.addSpectrum(s)

    # First load data and cache to disk
    pyopenms.MzMLFile().store("tfile_idx.mzML", mse)

    # Now load data
    ih = pyopenms.IndexedMzMLHandler("tfile_idx.mzML")

    assert ih.getNrChromatograms() ==0
    assert ih.getNrSpectra() == 1

    s = ih.getMSSpectrumById(0)
    s2 = ih.getSpectrumById(0)


@report
def testIDMapper():
    """
    @tests: IDMapper
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
    @tests: IdXMLFile
     IdXMLFile.__init__
     IdXMLFile.load
     IdXMLFile.store
    """
    assert pyopenms.IdXMLFile().load is not None
    assert pyopenms.IdXMLFile().store is not None

@report
def test_peptide_identifications_to_df():
    # convert to dataframe
    peps = []

    p = pyopenms.PeptideIdentification()
    p.setRT(1243.56)
    p.setMZ(440.0)
    p.setScoreType("ScoreType")
    p.setHigherScoreBetter(False)
    p.setIdentifier("IdentificationRun1")

    h = pyopenms.PeptideHit()
    h.setScore(1.0)
    h.setCharge(2)
    h.setMetaValue("StringMetaValue", "Value")
    h.setMetaValue("IntMetaValue", 2)
    e1 = pyopenms.PeptideEvidence()
    e1.setProteinAccession("sp|Accession1")
    e1.setStart(123)
    e1.setEnd(141)
    e2 = pyopenms.PeptideEvidence()
    e2.setProteinAccession("sp|Accession2")
    e2.setStart(12)
    e2.setEnd(24)
    h.setPeptideEvidences([e1, e2])
    p.insertHit(h)

    peps.append(p)

    p1 = pyopenms.PeptideIdentification()
    p1.setRT(1243.56)
    p1.setMZ(240.0)
    p1.setScoreType("ScoreType")
    p1.setHigherScoreBetter(False)
    p1.setIdentifier("IdentificationRun2")

    peps.append(p1)

    assert pyopenms.peptide_identifications_to_df(peps).shape == (2,12)
    assert pyopenms.peptide_identifications_to_df(peps, decode_ontology=False).shape == (2,12)
    assert pyopenms.peptide_identifications_to_df(peps)['protein_accession'][0] == 'sp|Accession1,sp|Accession2'
    assert pyopenms.peptide_identifications_to_df(peps, export_unidentified=False).shape == (1,12)

    # update from dataframe
    df = pyopenms.peptide_identifications_to_df(peps)
    df.loc[0, "ScoreType"] = 10.0
    peps = pyopenms.update_scores_from_df(peps, df, "ScoreType")
    assert peps[0].getHits()[0].getScore() == 10.0

@report
def testPepXMLFile():
    """
    @tests: PepXMLFile
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
    @tests: ProtXMLFile
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
    @tests: DTA2DFile
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
    @tests: DTAFile
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
    @tests: EDTAFile
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
    @tests: KroenikFile
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
    @tests: MSPFile
     MSPFile.__init__
    """
    f = pyopenms.MSPFile()

    # assert pyopenms.KroenikFile().load is not None
    # assert pyopenms.KroenikFile().store is not None

@report
def testMzIdentMLFile():
    """
    @tests: MzIdentMLFile
     MzIdentMLFile.__init__
    """
    f = pyopenms.MzIdentMLFile()

    assert pyopenms.MzIdentMLFile().load is not None
    assert pyopenms.MzIdentMLFile().store is not None
    assert pyopenms.MzIdentMLFile().isSemanticallyValid is not None


@report
def testMzTabFile():
    """
    @tests: MzTabFile
     MzTabFile.__init__
    """
    f = pyopenms.MzTabFile()

    # assert pyopenms.MzTabFile().store is not None

@report
def testMzTab():
    """
    @tests: MzTab
     MzTab.__init__
    """
    # f = pyopenms.MzTab()

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
    @tests: LogType
     LogType.CMD
     LogType.GUI
     LogType.NONE
     """
    assert isinstance(pyopenms.LogType.CMD, int)
    assert isinstance(pyopenms.LogType.GUI, int)
    assert isinstance(pyopenms.LogType.NONE, int)

# performance measurement helper for XIC and peak extraction
import time
import random
from typing import Tuple, List

def generate_random_ranges(exp: pyopenms.MSExperiment, 
                         n_ranges: int,
                         rt_width: float,
                         mz_width: float) -> np.ndarray:
    """
    Generate random ranges within the experiment bounds.
    
    Args:
        exp: MSExperiment object
        n_ranges: Number of ranges to generate
        rt_width: Width of RT window
        mz_width: Width of m/z window
        
    Returns:
        numpy array of shape (n_ranges, 4) with [mz_min, mz_max, rt_min, rt_max]
    """

    # Set the seed
    np.random.seed(4711)
    
    min_mz = exp.getMinMZ()
    max_mz = exp.getMaxMZ()
    min_rt = exp.getMinRT()
    max_rt = exp.getMaxRT()    
    print(f"spectra min_mz: {min_mz}, max_mz: {max_mz}, min_rt: {min_rt}, max_rt: {max_rt}")

    # Generate random centers
    mz_centers = np.random.uniform(min_mz + mz_width/2, max_mz - mz_width/2, n_ranges)
    rt_centers = np.random.uniform(min_rt + rt_width/2, max_rt - rt_width/2, n_ranges)
    
    # Create ranges array
    ranges = np.zeros((n_ranges, 4))
    ranges[:, 0] = mz_centers - mz_width/2  # mz_min
    ranges[:, 1] = mz_centers + mz_width/2  # mz_max
    ranges[:, 2] = rt_centers - rt_width/2  # rt_min
    ranges[:, 3] = rt_centers + rt_width/2  # rt_max
    
    return ranges

def extraction_performance_test(exp: pyopenms.MSExperiment, ms_level: int) -> Tuple[float, float]:
    """
    Run performance test for both aggregateFromMatrix and extractXICsFromMatrix
    
    Args:
        exp: MSExperiment object
        
    Returns:
        Tuple of (aggregate_time, xic_time)
    """
    # Generate 1_000,000 random ranges
    print("\nGenerating random ranges...")
    ranges = generate_random_ranges(exp, 
                                  n_ranges=1_000_000,
                                  rt_width=60.0,
                                  mz_width=0.01)
    
    # Convert to MatrixDouble
    ranges_matrix = pyopenms.MatrixDouble.fromNdArray(ranges)
    
    # Time aggregateFromMatrix
    print("Running aggregateFromMatrix...")
    start_time = time.time()
    _ = exp.aggregateFromMatrix(ranges_matrix, ms_level, b"sum")
    aggregate_time = time.time() - start_time
    
    # Time extractXICsFromMatrix
    print("Running extractXICsFromMatrix...")
    start_time = time.time()
    _ = exp.extractXICsFromMatrix(ranges_matrix, ms_level, b"sum")
    xic_time = time.time() - start_time
    
    return aggregate_time, xic_time

@report
def testMSExperiment():
    """
    @tests: MSExperiment
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
     MSExperiment.get2DPeakDataLong
     MSExperiment.get_df
     MSExperiment.get_massql_df
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
    _testStrOutput(mse.getLoadedFilePath())
    assert isinstance(mse.getMinIntensity(), float)
    assert isinstance(mse.getMaxIntensity(), float)

    mse.setLoadedFilePath("")
    assert mse.size() == 0

    mse.getIdentifier()
    mse.getLoadedFileType()
    mse.getLoadedFilePath()

    spec = pyopenms.MSSpectrum()
    data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
    data_i = np.array( [50.0, 80.0] ).astype(np.float32)
    spec.set_peaks( [data_mz,data_i] )

    mse.addSpectrum(spec)
    assert mse.size() == 1

    assert mse[0] is not None

    mse.updateRanges()
    rt, mz, inty = mse.get2DPeakDataLong(mse.getMinRT(), mse.getMaxRT(), mse.getMinMZ(), mse.getMaxMZ(), 1)
    assert rt.shape[0] == 2
    assert mz.shape[0] == 2
    assert inty.shape[0] == 2

    assert isinstance(list(mse), list)

    assert mse == mse
    assert not mse != mse

    assert mse.getSize() >= 0
    assert int(mse.isSorted()) in (0,1)

    mse2 = copy.copy(mse)

    assert mse.getSize() == mse2.getSize()
    assert mse2 == mse

    exp = pyopenms.MSExperiment()

    for i in range(5):
        s = pyopenms.MSSpectrum()
        s.setRT(i)
        s.setMSLevel(1 if i % 2 == 0 else 2)

        for mz in (500, 600):
            p = pyopenms.Peak1D()
            p.setMZ(mz + i)
            p.setIntensity(i + 10)
            s.push_back(p)

        exp.addSpectrum(s)

    assert exp.get_df().shape == (5, 4)
    assert exp.get_df(ms_levels=[1]).shape == (3, 4)
    assert exp.get_df(ms_levels=[2]).shape == (2, 4)

    assert exp.get_df(long=True).shape == (10, 4)
    assert exp.get_df(long = True, ms_levels=[1]).shape == (6, 4)
    assert exp.get_df(long=True, ms_levels=[2]).shape == (4, 4)

    pyopenms.MzMLFile().load(os.path.join(os.environ['OPENMS_DATA_PATH'], 'examples/FRACTIONS/BSA1_F1.mzML'), exp)

    ms1_df, ms2_df = exp.get_massql_df()
    assert ms1_df.shape == (140055, 7)
    assert np.allclose(ms2_df.head(), pd.read_csv(os.path.join(os.environ['OPENMS_DATA_PATH'], 'examples/FRACTIONS/BSA1_F1_MS2_MassQL.tsv'), sep='\t'))

    pyopenms.MzMLFile().load(os.path.join(os.environ['OPENMS_DATA_PATH'], 'examples/FRACTIONS/BSA1_F1_ION.mzML'), exp)
    df = exp.get_ion_df()
    assert np.allclose(df.head(), pd.read_csv(os.path.join(os.environ['OPENMS_DATA_PATH'], 'examples/FRACTIONS/BSA1_F1_MS1_ION.tsv'), sep='\t'))

    ms1_df, ms2_df = exp.get_massql_df(ion_mobility=True)
    assert ms1_df.shape == (332620, 8)
    assert np.allclose(ms1_df.head(), pd.read_csv(os.path.join(os.environ['OPENMS_DATA_PATH'], 'examples/FRACTIONS/BSA1_F1_MS1_MassQL_ION.tsv'), sep='\t'))

    #####################################################################################
    # test fast aggregation and XIC extraction using ranges
    pyopenms.MzMLFile().load(os.path.join(os.environ['OPENMS_DATA_PATH'], 'examples/FRACTIONS/BSA1_F1.mzML'), exp)    
    exp.updateRanges()

    ############################################################################
    # Uncomment to run performance tests
    # print("\nStarting performance tests MS1...")
    # aggregate_time, xic_time = extraction_performance_test(exp, 1) # Run performance tests on MS level 1   
    # print("\nPerformance Results:")
    # print(f"aggregateFromMatrix time for 1M ranges: {aggregate_time:.2f} seconds")
    # print(f"extractXICsFromMatrix time for 1M ranges: {xic_time:.2f} seconds")
    # print(f"Ratio (XIC/aggregate): {xic_time/aggregate_time:.2f}")

    # print("\nStarting performance tests MS2...")
    # aggregate_time, xic_time = extraction_performance_test(exp, 2) # Run performance tests on MS level 1   
    # print("\nPerformance Results:")
    # print(f"aggregateFromMatrix time for 1M ranges: {aggregate_time:.2f} seconds")
    # print(f"extractXICsFromMatrix time for 1M ranges: {xic_time:.2f} seconds")
    # print(f"Ratio (XIC/aggregate): {xic_time/aggregate_time:.2f}")
    # assert false # needed to output the results

    # eluting peptide feature at these coordinates
    rt_min = 1730.0
    rt_max = 1770.0
    iso1_mz = 443.711 # monoisotopic peak
    iso2_mz = 444.212 # first isotopic peak
    iso3_mz = 444.713 # second isotopic peak

    no_iso1_mz = 444.000 # no peaks in this area
    no_iso2_mz = 444.403 # some noise peaks in this area


    # Create ranges matrix with structure:
    # [[mz_min, mz_max, rt_min, rt_max], ...]
    ranges_matrix = pyopenms.MatrixDouble.fromNdArray(
         np.array([
        # Expected isotope peaks
        [iso1_mz - 0.1, iso1_mz + 0.1, rt_min, rt_max],
        [iso2_mz - 0.1, iso2_mz + 0.1, rt_min, rt_max],
        [iso3_mz - 0.1, iso3_mz + 0.1, rt_min, rt_max],
        # Regions where we don't expect peaks
        [no_iso1_mz - 0.1, no_iso1_mz + 0.1, rt_min, rt_max],
        [no_iso2_mz - 0.1, no_iso2_mz + 0.1, rt_min, rt_max]
    ])
    )
    
    # Print ranges matrix for verification
    print("Ranges Matrix:")
    print (ranges_matrix.get_matrix())

    agg_result = exp.aggregateFromMatrix(ranges_matrix, 1, b"sum")    
    agg_result_array = np.array(agg_result)    # Convert result to numpy array for easier testing

    print("\nAggregation Results:")
    print(agg_result)

    # Basic shape checks
    assert len(agg_result_array) == ranges_matrix.rows(), f"Expected {ranges_matrix.rows()} results, got {len(agg_result_array)}"
    
    # Check that regions with expected peaks have higher values
    isotope_intensities = agg_result_array[:3]  # First three rows are isotope peaks
    no_peak_intensities = agg_result_array[3:]  # Last two rows are regions without peaks
    
    print(isotope_intensities)
    print(no_peak_intensities)
 
    print("\nIsotope Peak Arrays and Their Sums:")
    isotope_sums = []
    for i, intensity_array in enumerate(isotope_intensities):
        array_sum = np.sum(intensity_array)
        isotope_sums.append(array_sum)
        print(f"Isotope {i+1} array: {intensity_array}")
        print(f"Isotope {i+1} sum: {array_sum}")
    
    print("\nNo-Peak Region Arrays and Their Sums:")
    no_peak_sums = []
    for i, intensity_array in enumerate(no_peak_intensities):
        array_sum = np.sum(intensity_array)
        no_peak_sums.append(array_sum)
        print(f"No-peak region {i+1} array: {intensity_array}")
        print(f"No-peak region {i+1} sum: {array_sum}")

    EXPECTED_ISO_SUMS = [24680058.50, 11043987.94, 3141677.76] 
    EXPECTED_NO_PEAK_SUMS = [0.0, 12322.06]  
    
    # Test each isotope array sum
    for i, (actual_sum, expected_sum) in enumerate(zip(isotope_sums, EXPECTED_ISO_SUMS)):
        assert np.isclose(actual_sum, expected_sum, rtol=1e-5), \
            f"Expected isotope {i+1} sum {expected_sum}, got {actual_sum}"
    
    # Test each no-peak array sum
    for i, (actual_sum, expected_sum) in enumerate(zip(no_peak_sums, EXPECTED_NO_PEAK_SUMS)):
        assert np.isclose(actual_sum, expected_sum, rtol=1e-5), \
            f"Expected no-peak region {i+1} sum {expected_sum}, got {actual_sum}"


    xic_result = exp.extractXICsFromMatrix(ranges_matrix, 1, b"sum")
    print("\nXIC Results:")
    for chrom in xic_result:
        print(f"m/z: {chrom.getProduct().getMZ()}")
        print(f"Size: {chrom.size()}")

    xic_details = []
    for i, chrom in enumerate(xic_result):
        # Get basic XIC information
        mz = chrom.getProduct().getMZ()
        size = chrom.size()
        
        # Get intensity values from the chromatogram
        intensities = [point.getIntensity() for point in chrom]
        rt_values = [point.getRT() for point in chrom]
        total_intensity = sum(intensities)
        
        details = {
            'mz': mz,
            'size': size,
            'total_intensity': total_intensity,
            'rt_range': (min(rt_values), max(rt_values)) if rt_values else (None, None)
        }
        xic_details.append(details)
        
        print(f"\nXIC {i+1}:")
        print(f"m/z: {mz}")
        print(f"Number of points: {size}")
        print(f"Total intensity: {total_intensity}")
        print(f"RT range: {details['rt_range']}")

    # XIC expectations
    EXPECTED_XIC_SIZES = [24, 24, 24, 24, 24]  # Expected number of points in each XIC
    EXPECTED_XIC_TOTAL_INTENSITIES = [24680058.50, 11043987.94, 3141677.76, 0.0, 12322.06]  # Expected total intensity for each XIC
    
    # Test XIC results
    for i, (details, exp_size, exp_total) in enumerate(zip(
            xic_details, EXPECTED_XIC_SIZES, EXPECTED_XIC_TOTAL_INTENSITIES)):
        
        assert details['size'] == exp_size, \
            f"XIC {i+1}: Expected {exp_size} points, got {details['size']}"
        
        assert np.isclose(details['total_intensity'], exp_total, rtol=1e-5), \
            f"XIC {i+1}: Expected total intensity {exp_total}, got {details['total_intensity']}"
                
        # Verify RT range falls within expected bounds
        if details['rt_range'][0] is not None:
            assert rt_min <= details['rt_range'][0] <= rt_max, \
                f"XIC {i+1}: Start RT {details['rt_range'][0]} outside expected range [{rt_min}, {rt_max}]"
            assert rt_min <= details['rt_range'][1] <= rt_max, \
                f"XIC {i+1}: End RT {details['rt_range'][1]} outside expected range [{rt_min}, {rt_max}]"


@report
def testMSSpectrum():
    """
    @tests: MSSpectrum
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

    assert isinstance(spec.getMinMZ(), float)
    assert isinstance(spec.getMaxMZ(), float)
    assert isinstance(spec.getMinIntensity(), float)
    assert isinstance(spec.getMaxIntensity(), float)

    assert spec == spec
    assert not spec != spec

    mz, ii = spec.get_peaks()
    assert len(mz) == len(ii)
    assert len(mz) == 1

    spec.set_peaks((mz, ii))
    mz0, ii0 = spec.get_peaks()
    assert mz0 == mz
    assert ii0 == ii

    assert int(spec.isSorted()) in (0,1)

    spec.clear(False)
    p = pyopenms.Peak1D()
    p.setMZ(1000.0)
    p.setIntensity(200.0)
    spec.push_back(p)
    p = pyopenms.Peak1D()
    p.setMZ(2000.0)
    p.setIntensity(400.0)
    spec.push_back(p)

    mz, ii = spec.get_peaks()
    assert spec[0].getMZ() == 1000.0
    assert spec[1].getMZ() == 2000.0
    assert spec[0].getIntensity() == 200.0
    assert spec[1].getIntensity() == 400.0
    assert mz[0] == 1000.0
    assert mz[1] == 2000.0
    assert ii[0] == 200.0
    assert ii[1] == 400.0

    spec.setMSLevel(2)
    spec.setNativeID('scan=1')
    prec = pyopenms.Precursor()
    prec.setMZ(100.0)
    prec.setCharge(1)
    spec.setPrecursors([prec])
    spec.setMetaValue('total ion current', 600)
    pepid = pyopenms.PeptideIdentification()
    hit = pyopenms.PeptideHit(1.0, 1, 0, pyopenms.AASequence.fromString('A'))
    pepid.setHits([hit])
    spec.setPeptideIdentifications([pepid])

    data = np.array( [5, 8] ).astype(np.float32)
    f_da = [ pyopenms.FloatDataArray() ]
    f_da[0].set_data(data)
    f_da[0].setName("Ion Mobility")
    spec.setFloatDataArrays( f_da )
    spec.setDriftTimeUnit( pyopenms.DriftTimeUnit.MILLISECOND )

    s_da = pyopenms.StringDataArray()
    for s in ['b3+', 'y4+']:
        s_da.push_back(s)
    s_da.setName("IonNames")
    spec.setStringDataArrays([s_da])

    df = spec.get_df()
    assert df.shape == (2, 11)
    assert df.loc[0, 'mz'] == 1000.0
    assert df.loc[0, 'intensity'] == 200.0
    assert df.loc[0, 'ion_mobility'] == 5.0
    assert df.loc[0, 'ion_mobility_unit'] == 'ms'
    assert df.loc[0, 'ms_level'] == 2
    assert df.loc[0, 'precursor_mz'] == 100.0
    assert df.loc[0, 'precursor_charge'] == 1
    assert df.loc[0, 'native_id'] == 'scan=1'
    assert df.loc[0, 'ion_annotation'] == 'b3+'
    assert df.loc[0, 'sequence'] == 'A'
    assert df.loc[0, 'total ion current'] == 600

    spec.clear(False)
    data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
    data_i = np.array( [50.0, 80.0] ).astype(np.float32)
    spec.set_peaks( [data_mz,data_i] )

    mz, ii = spec.get_peaks()
    assert spec[0].getMZ() == 5.0
    assert spec[1].getMZ() == 8.0
    assert spec[0].getIntensity() == 50.0
    assert spec[1].getIntensity() == 80.0
    assert mz[0] == 5.0
    assert mz[1] == 8.0
    assert ii[0] == 50.0
    assert ii[1] == 80.0

    # Fast
    spec.clear(False)
    data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
    data_i = np.array( [50.0, 80.0] ).astype(np.float64)
    spec.set_peaks( [data_mz,data_i] )

    mz, ii = spec.get_peaks()
    assert spec[0].getMZ() == 5.0
    assert spec[1].getMZ() == 8.0
    assert spec[0].getIntensity() == 50.0
    assert spec[1].getIntensity() == 80.0
    assert mz[0] == 5.0
    assert mz[1] == 8.0
    assert ii[0] == 50.0
    assert ii[1] == 80.0

    # Slow
    spec.clear(False)
    data_mz = np.array( [5.0, 8.0] ).astype(np.float32)
    data_i = np.array( [50.0, 80.0] ).astype(np.float32)
    spec.set_peaks( [data_mz,data_i] )

    mz, ii = spec.get_peaks()
    assert spec[0].getMZ() == 5.0
    assert spec[1].getMZ() == 8.0
    assert spec[0].getIntensity() == 50.0
    assert spec[1].getIntensity() == 80.0
    assert mz[0] == 5.0
    assert mz[1] == 8.0
    assert ii[0] == 50.0
    assert ii[1] == 80.0

    ###################################
    # get data arrays
    ###################################
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


    spec = pyopenms.MSSpectrum()
    assert len(spec.getIntegerDataArrays()) == 0
    int_da = [ pyopenms.IntegerDataArray() ]
    int_da[0].push_back(5)
    int_da[0].push_back(6)
    int_da.append( pyopenms.IntegerDataArray() )
    int_da[1].push_back(8)
    spec.setIntegerDataArrays( int_da )
    assert len(spec.getIntegerDataArrays()) == 2
    assert spec.getIntegerDataArrays()[0][0] == 5
    assert spec.getIntegerDataArrays()[1][0] == 8

    spec = pyopenms.MSSpectrum()
    data = np.array( [5, 8, 42] ).astype(np.intc)
    int_da = [ pyopenms.IntegerDataArray() ]
    int_da[0].set_data(data)
    spec.setIntegerDataArrays( int_da )
    assert len(spec.getIntegerDataArrays()) == 1
    assert spec.getIntegerDataArrays()[0][0] == 5
    assert spec.getIntegerDataArrays()[0][2] == 42
    assert len(int_da[0].get_data() ) == 3

    spec = pyopenms.MSSpectrum()
    assert len(spec.getFloatDataArrays()) == 0
    f_da = [ pyopenms.FloatDataArray() ]
    f_da[0].push_back(5.0)
    f_da[0].push_back(6.0)
    f_da.append( pyopenms.FloatDataArray() )
    f_da[1].push_back(8.0)
    spec.setFloatDataArrays( f_da )
    assert len(spec.getFloatDataArrays()) == 2.0
    assert spec.getFloatDataArrays()[0][0] == 5.0
    assert spec.getFloatDataArrays()[1][0] == 8.0

    spec = pyopenms.MSSpectrum()
    data = np.array( [5, 8, 42] ).astype(np.float32)
    f_da = [ pyopenms.FloatDataArray() ]
    f_da[0].set_data(data)
    spec.setFloatDataArrays( f_da )
    assert len(spec.getFloatDataArrays()) == 1
    assert spec.getFloatDataArrays()[0][0] == 5.0
    assert spec.getFloatDataArrays()[0][2] == 42.0
    assert len(f_da[0].get_data() ) == 3

    spec = pyopenms.MSSpectrum()
    dfunit = spec.getDriftTimeUnit()
    assert pyopenms.DriftTimeUnit().getMapping()[dfunit]  == "NONE"
    assert dfunit == pyopenms.DriftTimeUnit.NONE
    assert spec.getDriftTimeUnitAsString() == '<NONE>'
    spec.setDriftTimeUnit( pyopenms.DriftTimeUnit.MILLISECOND )

    dfunit = spec.getDriftTimeUnit()
    assert dfunit == pyopenms.DriftTimeUnit.MILLISECOND
    assert pyopenms.DriftTimeUnit().getMapping()[dfunit]  == "MILLISECOND"
    assert spec.getDriftTimeUnitAsString() == 'ms'

    spec = pyopenms.MSSpectrum()
    spec.setDriftTime(6.0)
    assert spec.getDriftTime() == 6.0

    spec = pyopenms.MSSpectrum()
    assert not spec.containsIMData()
    data = np.array( [5, 8, 42] ).astype(np.float32)
    f_da = [ pyopenms.FloatDataArray() ]
    f_da[0].set_data(data)
    f_da[0].setName("Ion Mobility")
    spec.setFloatDataArrays( f_da )
    assert spec.containsIMData()
    assert spec.getIMData()[0] == 0
    f_da = [ pyopenms.FloatDataArray(), pyopenms.FloatDataArray() ]
    f_da[0].setName("test")
    f_da[0].set_data(data)
    f_da[1].set_data(data)
    f_da[1].setName("Ion Mobility")
    spec.setFloatDataArrays( f_da )
    assert spec.containsIMData()
    assert spec.getIMData()[0] == 1

    # Ensure that "set_peaks()" doesnt clear the float data arrays
    spec = pyopenms.MSSpectrum()
    data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
    data_i = np.array( [50.0, 80.0] ).astype(np.float32)
    f_da = [ pyopenms.FloatDataArray() ]
    f_da[0].set_data(data)
    f_da[0].setName("Ion Mobility")
    spec.setFloatDataArrays( f_da )
    spec.set_peaks( [data_mz,data_i] )
    assert spec.containsIMData()
    assert spec.getIMData()[0] == 0
    assert len(spec.getFloatDataArrays()) == 1

    f = spec.getFloatDataArrays()[0]
    assert len(f.get_data()) == 3
    assert f.get_data()[0] == 5
    assert spec.size() == len(data_mz)
    assert spec.size() == len(data_i)

@report
def testStringDataArray():
    """
    @tests: StringDataArray
     """
    da = pyopenms.StringDataArray()
    assert da.size() == 0
    da.push_back("hello")
    da.push_back("world")
    assert da.size() == 2
    assert da[0] == b"hello"
    assert da[1] == b"world"
    da[1] = b"hello world"
    assert da[1] == b"hello world", da[1]
    da.clear()
    assert da.size() == 0
    da.push_back("hello")
    assert da.size() == 1
    da.resize(3)
    da[0] = b"hello"
    da[1] = b""
    da[2] = b"world"
    assert da.size() == 3

@report
def testIntegerDataArray():
    """
    @tests: IntegerDataArray
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

    q = da.get_data()
    q = np.append(q, 4).astype(np.intc)
    da.set_data(q)
    assert da.size() == 4

@report
def testFloatDataArray():
    """
    @tests: FloatDataArray
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

    q = da.get_data()
    q = np.append(q, 4.0).astype(np.float32)
    da.set_data(q)
    assert da.size() == 4

@report
def testMSChromatogram():
    """
    @tests: MSChromatogram
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

    chrom.setName("chrom")
    assert chrom.getName() == "chrom"

    p = pyopenms.ChromatogramPeak()
    p.setRT(1000.0)
    p.setIntensity(200.0)

    chrom.push_back(p)
    assert chrom.size() == 1
    assert chrom[0] == p

    chrom.updateRanges()
    assert isinstance(chrom.findNearest(0.0), int)

    assert isinstance(chrom.getMinRT(), float)
    assert isinstance(chrom.getMaxRT(), float)
    assert isinstance(chrom.getMinIntensity(), float)
    assert isinstance(chrom.getMaxIntensity(), float)

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

    chrom.clear(False)
    p = pyopenms.ChromatogramPeak()
    p.setRT(1000.0)
    p.setIntensity(200.0)
    chrom.push_back(p)
    p = pyopenms.ChromatogramPeak()
    p.setRT(2000.0)
    p.setIntensity(400.0)
    chrom.push_back(p)

    mz, ii = chrom.get_peaks()
    assert chrom[0].getRT() == 1000.0
    assert chrom[1].getRT() == 2000.0
    assert chrom[0].getIntensity() == 200.0
    assert chrom[1].getIntensity() == 400.0
    assert mz[0] == 1000.0
    assert mz[1] == 2000.0
    assert ii[0] == 200.0
    assert ii[1] == 400.0

    chrom.setNativeID('chrom_0')
    chrom.setMetaValue('FWHM', 5.0)
    prec = pyopenms.Precursor()
    prec.setMZ(100.0)
    prec.setCharge(1)
    chrom.setPrecursor(prec)
    prod = pyopenms.Product()
    prod.setMZ(50.0)
    chrom.setProduct(prod)
    chrom.setComment('comment')

    df = chrom.get_df()
    assert df.shape == (2, 9)
    assert df.loc[0, 'time'] == 1000.0
    assert df.loc[1, 'intensity'] == 400
    assert df.loc[0, 'chromatogram_type'] == 'MASS_CHROMATOGRAM'
    assert df.loc[1, 'precursor_mz'] == 100.0
    assert df.loc[0, 'precursor_charge'] == 1
    assert df.loc[1, 'product_mz'] == 50.0
    assert df.loc[0, 'comment'] == 'comment'
    assert df.loc[1, 'native_id'] == 'chrom_0'
    assert df.loc[0, 'FWHM'] == 5

    chrom.clear(False)
    data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
    data_i = np.array( [50.0, 80.0] ).astype(np.float32)
    chrom.set_peaks( [data_mz,data_i] )

    mz, ii = chrom.get_peaks()
    assert chrom[0].getRT() == 5.0
    assert chrom[1].getRT() == 8.0
    assert chrom[0].getIntensity() == 50.0
    assert chrom[1].getIntensity() == 80.0
    assert mz[0] == 5.0
    assert mz[1] == 8.0
    assert ii[0] == 50.0
    assert ii[1] == 80.0

    # Fast
    chrom.clear(False)
    data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
    data_i = np.array( [50.0, 80.0] ).astype(np.float64)
    chrom.set_peaks( [data_mz,data_i] )

    mz, ii = chrom.get_peaks()
    assert chrom[0].getRT() == 5.0
    assert chrom[1].getRT() == 8.0
    assert chrom[0].getIntensity() == 50.0
    assert chrom[1].getIntensity() == 80.0
    assert mz[0] == 5.0
    assert mz[1] == 8.0
    assert ii[0] == 50.0
    assert ii[1] == 80.0

    # Slow
    chrom.clear(False)
    data_mz = np.array( [5.0, 8.0] ).astype(np.float32)
    data_i = np.array( [50.0, 80.0] ).astype(np.float32)
    chrom.set_peaks( [data_mz,data_i] )

    mz, ii = chrom.get_peaks()
    assert chrom[0].getRT() == 5.0
    assert chrom[1].getRT() == 8.0
    assert chrom[0].getIntensity() == 50.0
    assert chrom[1].getIntensity() == 80.0
    assert mz[0] == 5.0
    assert mz[1] == 8.0
    assert ii[0] == 50.0
    assert ii[1] == 80.0

    # test float data
    chrom = pyopenms.MSChromatogram()
    data = np.array( [5, 8, 42] ).astype(np.float32)
    f_da = [ pyopenms.FloatDataArray() ]
    f_da[0].set_data(data)
    f_da[0].setName("Test Data")
    chrom.setFloatDataArrays( f_da )
    assert len(chrom.getFloatDataArrays()) == 1
    f = chrom.getFloatDataArrays()[0]
    assert f.get_data()[0] == 5
    assert f.getName() == "Test Data"

    # Ensure that "set_peaks()" doesnt clear the float data arrays
    chrom = pyopenms.MSChromatogram()
    chrom.setFloatDataArrays( f_da )
    chrom.set_peaks( [data_mz,data_i] )
    assert len(chrom.getFloatDataArrays()) == 1

    f = chrom.getFloatDataArrays()[0]
    assert len(f.get_data()) == 3
    assert f.get_data()[0] == 5
    assert chrom.size() == len(data_mz)
    assert chrom.size() == len(data_i)


@report
def testMRMFeature():
    """
    @tests: MRMFeature
      MRMFeature.__init__
      MRMFeature.addScore
      MRMFeature.getScore
     """
    mrmfeature = pyopenms.MRMFeature()
    f = pyopenms.Feature()

    fs = mrmfeature.getFeatures()
    assert len(fs) == 0

    mrmfeature.addFeature(f, "myFeature")
    fs = mrmfeature.getFeatures()
    assert len(fs) == 1
    assert mrmfeature.getFeature("myFeature") is not None
    slist = []
    mrmfeature.getFeatureIDs(slist)
    assert len(slist) == 1

    mrmfeature.addPrecursorFeature(f, "myFeature_Pr0")
    slist = []
    mrmfeature.getPrecursorFeatureIDs(slist)
    assert len(slist) == 1
    assert mrmfeature.getPrecursorFeature("myFeature_Pr0") is not None

    s = mrmfeature.getScores()
    assert abs(s.yseries_score - 0.0) < 1e-4
    s.yseries_score = 4.0
    mrmfeature.setScores(s)
    s2 = mrmfeature.getScores()
    assert abs(s2.yseries_score - 4.0) < 1e-4

@report
def testConfidenceScoring():
    """
    @tests: ConfidenceScoring
      ConfidenceScoring.__init__
     """
    scoring = pyopenms.ConfidenceScoring()

@report
def testMRMDecoy():
    """
    @tests: MRMDecoy
      MRMDecoy.__init__
     """
    mrmdecoy = pyopenms.MRMDecoy()
    assert mrmdecoy is not None

    assert pyopenms.MRMDecoy().generateDecoys is not None

@report
def testMRMTransitionGroup():
    """
    @tests: MRMTransitionGroup
    """
    mrmgroup = pyopenms.MRMTransitionGroupCP()
    assert mrmgroup is not None

    mrmgroup.setTransitionGroupID("this_id")
    assert mrmgroup.getTransitionGroupID() == "this_id"

    assert len(mrmgroup.getTransitions()) == 0
    mrmgroup.addTransition(pyopenms.ReactionMonitoringTransition(), "tr1")
    assert len(mrmgroup.getTransitions()) == 1

    # add data for testing df output
    ## test chromatogram df
    rt, intensity = [[1.0], [5]]
    chrom = pyopenms.MSChromatogram()
    chrom.set_peaks([rt, intensity])
    chrom.setNativeID("tr1")
    mrmgroup.addChromatogram(chrom, 'tr1')

    df = mrmgroup.get_chromatogram_df()
    assert df.shape == (1, 8)
    assert df.loc[0, 'time'] == 1.0
    assert df.loc[0, 'intensity'] == 5
    assert df.loc[0, 'chromatogram_type'] == 'MASS_CHROMATOGRAM'
    assert df.loc[0, 'native_id'] == 'tr1'

    ## feature 1
    f1 = pyopenms.MRMFeature()
    f1.setRT(1.0)
    f1.setMetaValue(b'leftWidth', 0.5)
    f1.setMetaValue(b'rightWidth', 1.5)
    f1.setMetaValue(b'peak_apices_sum', 10.0)
    f1.setOverallQuality(0.5)
    f1.setUniqueId(1)
    f1.setIntensity(20.0)
    mrmgroup.addFeature(f1)

    # feature 2
    f2 = pyopenms.MRMFeature()
    f2.setRT(2.0)
    f2.setMetaValue(b'peak_apices_sum', 20.0)
    f2.setOverallQuality(1.0)
    f2.setUniqueId(2)
    f2.setIntensity(40.0)
    mrmgroup.addFeature(f2)

    df = mrmgroup.get_feature_df(meta_values=[b'leftWidth', b'rightWidth', b'peak_apices_sum'])
    assert df.shape == (2, 6)
    assert df.loc[1, 'leftWidth'] == 0.5
    assert df.loc[1, 'rightWidth'] == 1.5
    assert df.loc[1, 'peak_apices_sum'] == 10.0
    assert df.loc[1, 'intensity'] == 20.0
    assert df.loc[1, 'quality'] == 0.5
    assert df.loc[1, 'RT'] == 1.0

    assert np.isnan(df.loc[2, 'leftWidth'])
    assert np.isnan(df.loc[2, 'rightWidth'])
    assert df.loc[2, 'peak_apices_sum'] == 20.0
    assert df.loc[2, 'intensity'] == 40.0
    assert df.loc[2, 'quality'] == 1.0
    assert df.loc[2, 'RT'] == 2.0

    # If get "all" meta values should get the same result
    df = mrmgroup.get_feature_df(meta_values='all')
    assert df.shape == (2, 6)
    assert df.loc[1, 'leftWidth'] == 0.5
    assert df.loc[1, 'rightWidth'] == 1.5
    assert df.loc[1, 'peak_apices_sum'] == 10.0
    assert df.loc[1, 'intensity'] == 20.0
    assert df.loc[1, 'quality'] == 0.5
    assert df.loc[1, 'RT'] == 1.0

    assert np.isnan(df.loc[2, 'leftWidth'])
    assert np.isnan(df.loc[2, 'rightWidth'])
    assert df.loc[2, 'peak_apices_sum'] == 20.0
    assert df.loc[2, 'intensity'] == 40.0
    assert df.loc[2, 'quality'] == 1.0
    assert df.loc[2, 'RT'] == 2.0

@report
def testReactionMonitoringTransition():
    """
    @tests: ReactionMonitoringTransition
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
    @tests: MapAlignmentAlgorithmPoseClustering
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
    _testStrOutput(ma.getName())

    ma.setName(ma.getName())

    ma.getDefaults()
    ma.getParameters()

    ma.setParameters(ma.getDefaults())

    ma.setReference
    ma.align

    pyopenms.MapAlignmentTransformer.transformRetentionTimes

@report
def testMatrixDouble():

    """
    @tests: MatrixDouble
     MapAlignmentAlgorithmIdentification.__init__
     """

    m = pyopenms.MatrixDouble(3, 2, 0.0)
    for i in range(3):
        for j in range(2):
            m.setValue(i, j, i * 10.0 + j) 
    print(m)

    mv = m.get_matrix_as_view()
    print(mv)

    mc = m.get_matrix()
    print(mc)

    mat = m.get_matrix_as_view()

    N = 90
    m = pyopenms.MatrixDouble(N-1, N+2, 5.0)

    assert m.rows() == 89
    assert m.cols() == 92

    rows = N-1
    cols = N+2
    test = []
    for i in range(int(rows)):
        for j in range(int(cols)):
            test.append( m.getValue(i,j) )

    testm = np.asarray(test)
    testm = testm.reshape(rows, cols)

    assert sum(sum(testm)) == 40940.0
    assert sum(sum(testm)) == (N-1)*(N+2)*5

    matrix = m.get_matrix()
    assert sum(sum(matrix)) == 40940.0
    assert sum(sum(matrix)) == (N-1)*(N+2)*5

    matrix_view = m.get_matrix_as_view()
    assert sum(sum(matrix_view)) == 40940.0
    assert sum(sum(matrix_view)) == (N-1)*(N+2)*5


    # Column = 1 / Row = 2
    ## Now change a value:

    assert m.getValue(1, 2) == 5.0
    m.setValue(1, 2, 8.0)
    assert m.getValue(1, 2) == 8.0

    print(m)
    mat = m.get_matrix_as_view()
    print(mat)
    assert mat[1, 2] == 8.0

    mat = m.get_matrix()
    assert m.getValue(1, 2) == 8.0
    assert mat[1, 2] == 8.0

    # Whatever we change here gets changed in the raw data as well
    matrix_view = m.get_matrix_as_view()
    matrix_view[1, 6] = 11.0
    assert m.getValue(1, 6) == 11.0
    assert matrix_view[1, 6] == 11.0

    m = pyopenms.MatrixDouble()
    assert m.rows() == 0
    assert m.cols() == 0

    mat[3, 6] = 9.0
    m.set_matrix(mat)
    assert m.getValue(1, 2) == 8.0
    assert m.getValue(3, 6) == 9.0


@report
def testMapAlignmentIdentification():

    """
    @tests: MapAlignmentAlgorithmIdentification
     MapAlignmentAlgorithmIdentification.__init__
     """
    ma = pyopenms.MapAlignmentAlgorithmIdentification()

    assert pyopenms.MapAlignmentAlgorithmIdentification().align is not None
    assert pyopenms.MapAlignmentAlgorithmIdentification().setReference is not None

@report
def testMapAlignmentTransformer():

    """
    @tests: MapAlignmentTransformer
     MapAlignmentTransformer.__init__
     """
    ma = pyopenms.MapAlignmentTransformer()

    assert pyopenms.MapAlignmentTransformer().transformRetentionTimes is not None

@report
def testMxxxFile():
    """
    @tests: MzDataFile
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

    """
    mse = pyopenms.MSExperiment()
    s = pyopenms.MSSpectrum()
    mse.addSpectrum(s)

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

    myStr = pyopenms.String()
    fh.storeBuffer(myStr, mse)
    assert len(myStr.toString()) == 5269
    mse2 = pyopenms.MSExperiment()
    fh.loadBuffer(bytes(myStr), mse2)
    assert mse2 == mse
    assert mse2.size() == 1

    fh = pyopenms.MzXMLFile()
    _testProgressLogger(fh)
    fh.store("test.mzXML", mse)
    fh.load("test.mzXML", mse)
    fh.setOptions(fh.getOptions())



@report
def testParamXMLFile():

    """
    @tests: ParamXMLFile
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
    @tests: Peak1D
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
def testNumpressCoder():
    """
    """

    np = pyopenms.MSNumpressCoder()

    nc = pyopenms.NumpressConfig()
    nc.np_compression = np.NumpressCompression.LINEAR
    nc.estimate_fixed_point = True
    tmp = pyopenms.String()
    out = []
    inp =  [1.0, 2.0, 3.0]
    np.encodeNP(inp, tmp, True, nc)

    res = tmp.toString()
    assert len(res) != 0, len(res)
    assert res != "", res
    np.decodeNP(res, out, True, nc)
    assert len(out) == 3, (out, res)
    assert out == inp, out

    # Now try to use a simple Python string as input -> this will fail as we
    # cannot pass this by reference in C++
    res = ""
    try:
        np.encodeNP(inp, res, True, nc)
        has_error = False
    except AssertionError:
        has_error = True

    assert has_error

@report
def testNumpressConfig():
    """
    """

    n = pyopenms.MSNumpressCoder()
    np = pyopenms.NumpressConfig()
    np.np_compression = n.NumpressCompression.LINEAR
    assert np.np_compression == n.NumpressCompression.LINEAR
    np.numpressFixedPoint = 4.2
    np.numpressErrorTolerance = 4.2
    np.estimate_fixed_point = True
    np.linear_fp_mass_acc = 4.2
    np.setCompression("linear")

@report
def testBase64():
    """
    """

    b = pyopenms.Base64()
    out = pyopenms.String()
    inp =  [1.0, 2.0, 3.0]
    b.encode64(inp, b.ByteOrder.BYTEORDER_LITTLEENDIAN, out, False)
    res = out.toString()
    assert len(res) != 0
    assert res != ""

    convBack = []
    b.decode64(res, b.ByteOrder.BYTEORDER_LITTLEENDIAN, convBack, False)
    assert convBack == inp, convBack

    # For 32 bit
    out = pyopenms.String()
    b.encode32(inp, b.ByteOrder.BYTEORDER_LITTLEENDIAN, out, False)
    res = out.toString()
    assert len(res) != 0
    assert res != ""

    convBack = []
    b.decode32(res, b.ByteOrder.BYTEORDER_LITTLEENDIAN, convBack, False)
    assert convBack == inp, convBack

@report
def testPeakFileOptions():
    """
    @tests: PeakFileOptions
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
    @tests: MRMMapping
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
def testPeakPickerChromatogram():
    """
    @tests: PeakPickerChromatogram
     PeakPickerChromatogram.__init__
     PeakPickerChromatogram.pickChromatogram
    """

    p = pyopenms.PeakPickerChromatogram()
    assert p.pickChromatogram is not None

@report
def testPeakPickerHiRes():
    """
    @tests: PeakPickerHiRes
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
    @tests: PeakTypeEstimator
     PeakTypeEstimator.__init__
     PeakTypeEstimator.estimateType
    """

    pyopenms.PeakTypeEstimator().estimateType(pyopenms.MSSpectrum())

@report
def testPeptideHit():
    """
    @tests: PeptideHit
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

    ph = pyopenms.PeptideHit(1.0, 1, 0, pyopenms.AASequence.fromString("A"))
    _testMetaInfoInterface(ph)

    assert len(ph.getPeptideEvidences()) == 0
    assert ph.getPeptideEvidences() == []

    pe = pyopenms.PeptideEvidence()
    pe.setProteinAccession('B_id')

    ph.addPeptideEvidence(pe)
    assert len(ph.getPeptideEvidences()) == 1
    assert ph.getPeptideEvidences()[0].getProteinAccession() == 'B_id'

    ph.setPeptideEvidences([pe,pe])
    assert len(ph.getPeptideEvidences()) == 2
    assert ph.getPeptideEvidences()[0].getProteinAccession() == 'B_id'

    assert ph.getScore() == 1.0
    assert ph.getRank() == 1
    assert ph.getSequence().toString() == "A"

    ph.setScore(2.0)
    assert ph.getScore() == 2.0
    ph.setRank(30)
    assert ph.getRank() == 30
    ph.setSequence(pyopenms.AASequence.fromString("AAA"))
    assert ph.getSequence().toString() == "AAA"

    assert ph == ph
    assert not ph != ph

@report
def testPeptideEvidence():
    """
    @tests: PeptideEvidence
     PeptideEvidence.__init__
    """
    pe = pyopenms.PeptideEvidence()
    assert pe == pe
    assert not pe != pe

    pe.setProteinAccession('B_id')
    assert pe.getProteinAccession() == "B_id"

    pe.setAABefore(b'A')
    assert pe.getAABefore() == 'A'
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
    @tests: PeptideIdentification
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
    pe.setProteinAccession('B_id')

    ph = pyopenms.PeptideHit(1.0, 1, 0, pyopenms.AASequence.fromString("A"))
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
    _testStrOutput(pi.getScoreType())
    pi.setScoreType("A")
    assert isinstance(pi.isHigherScoreBetter(), int)
    _testStrOutput(pi.getIdentifier())
    pi.setIdentifier("id")
    pi.assignRanks()
    pi.sort()
    assert not pi.empty()

    pi.setSignificanceThreshold(6.0)


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
def testProteinHit():
    """
    @tests: ProteinHit
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
    @tests: ProteinIdentification
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
    @tests: RichPeak1D
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

@report
def testSpectrumSetting(s=pyopenms.SpectrumSettings()):
    """
    @tests: SpectrumSettings
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
                               pyopenms.SpectrumSettings.SpectrumType.CENTROID,
                               pyopenms.SpectrumSettings.SpectrumType.PROFILE]

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
    @tests: TransformationDescription
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
    @tests: TransformationModelInterpolated
     TransformationModelInterpolated.getDefaultParameters
     TransformationModelInterpolated.getParameters
     TransformationModelLinear.getDefaultParameters
     TransformationModelLinear.getParameters
     TransformationModelBSpline.getDefaultParameters
     TransformationModelBSpline.getParameters
     TransformationModelLowess.getDefaultParameters
     TransformationModelLowess.getParameters
     NB: THIS TEST STOPS AFTER THE FIRST FAILURE
    """
    for clz in [pyopenms.TransformationModelLinear,
                pyopenms.TransformationModelBSpline,
                pyopenms.TransformationModelInterpolated,
                pyopenms.TransformationModelLowess]:
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
    @tests: TransformationXMLFile
     TransformationXMLFile.__init__
     TransformationXMLFile.load
     TransformationXMLFile.store
    """
    fh = pyopenms.TransformationXMLFile()
    td = pyopenms.TransformationDescription()
    fh.store("test.transformationXML", td)
    fh.load("test.transformationXML", td, True)
    assert td.getDataPoints() == []

@report
def testIBSpectraFile():
    """
    @tests: IBSpectraFile
     IBSpectraFile.__init__
     IBSpectraFile.store
    """
    fh = pyopenms.IBSpectraFile()
    cmap = pyopenms.ConsensusMap()
    correctError = False
    try:
        fh.store( pyopenms.String("test.ibspectra.file"), cmap)
        assert False
    except RuntimeError:
        correctError = True

    assert correctError 

@report
def testSwathFile():
    """
    @tests: SwathFile
     SwathFile.__init__
     SwathFile.store
    """
    fh = pyopenms.SwathFile()

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
def testInspectInfile():
    """
    @tests: InspectInfile
     InspectInfile.__init__
    """
    inst = pyopenms.InspectInfile()

    assert inst.getModifications is not None
    mods = inst.getModifications()
    assert len(mods) == 0

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
        pos = decoder.findIndexListOffset("abcde", 100)
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

    regex_description = ""
    psi_id = ""
    xtandem_id = ""
    comet_id = 0
    omssa_id = 0
    e = pyopenms.DigestionEnzymeProtein("testEnzyme", "K", set([]), regex_description,
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
def testProteaseDB():
    edb = pyopenms.ProteaseDB()

    f = pyopenms.EmpiricalFormula()
    synonyms = set(["dummy", "other"])

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
    edb.hasElement(pyopenms.String("O"))

    e = edb.getElement(16)

    assert e.getName() == "Sulfur"
    assert e.getSymbol() == "S"
    assert e.getIsotopeDistribution()

    e2 = edb.getElement(pyopenms.String("O"))

    assert e2.getName() == "Oxygen"
    assert e2.getSymbol() == "O"
    assert e2.getIsotopeDistribution()

    # assume we discovered a new element
    e2 = edb.addElement(b"NewElement", b"NE", 300, {400 : 1.0}, {400 : 400.1}, False)
    e2 = edb.getElement(pyopenms.String("NE"))
    assert e2.getName() == "NewElement"

    # changing existing elements in tests might have side effects so we define a new element
    # add first new element
    e2 = edb.addElement(b"Kryptonite", b"@", 500, {999 : 0.7, 1000 : 0.3}, {999 : 999.01, 1000 : 1000.01}, False)
    e2 = edb.getElement(pyopenms.String("@"))
    assert e2.getName() == "Kryptonite"
    assert e2.getIsotopeDistribution()
    assert len(e2.getIsotopeDistribution().getContainer()) == 2
    assert abs(e2.getIsotopeDistribution().getContainer()[1].getIntensity() - 0.3) < 1e-5
    # replace element
    e2 = edb.addElement(b"Kryptonite", b"@", 500, {9999 : 1.0}, {9999 : 9999.1}, True)
    e2 = edb.getElement(pyopenms.String("@"))
    assert e2.getName() == "Kryptonite"
    assert e2.getIsotopeDistribution()
    assert len(e2.getIsotopeDistribution().getContainer()) == 1
    assert abs(e2.getIsotopeDistribution().getContainer()[0].getIntensity() - 1.0) < 1e-5
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

    assert rdb.hasResidue(s("Glycine"))
    glycine = rdb.getResidue(s("Glycine"))

    nrr = rdb.getNumberOfResidues()

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
    assert list(mods)[0].getFullId() == "Acetyl (N-term)"

    m = mdb.getModification(s("Carboxymethyl (C)"), "", pyopenms.ResidueModification.TermSpecificity.NUMBER_OF_TERM_SPECIFICITY)
    assert m.getFullId() == "Carboxymethyl (C)"

    m = mdb.getModification( s("Phosphorylation"), s("S"), pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m.getId() == "Phospho"

    # get out all mods (there should be many, some known ones as well!)
    mods = []
    m = mdb.getAllSearchModifications(mods)
    assert len(mods) > 100

    assert b"Phospho (S)" in mods
    assert b"Sulfo (S)" in mods
    assert not (b"Phospho" in mods)

    # search for specific modifications by mass
    m = mdb.getBestModificationByDiffMonoMass( 80.0, 1.0, "T", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getId() == "Phospho"
    assert m.getFullName() == "Phosphorylation"
    assert m.getUniModAccession() == "UniMod:21"

    m = mdb.getBestModificationByDiffMonoMass(80, 100, "T", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getId() == "Phospho"
    assert m.getFullName() == "Phosphorylation"
    assert m.getUniModAccession() == "UniMod:21"

    m = mdb.getBestModificationByDiffMonoMass(16, 1.0, "M", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    assert m is not None
    assert m.getId() == "Oxidation", m.getId()
    assert m.getFullName() == "Oxidation or Hydroxylation", m.getFullName()
    assert m.getUniModAccession() == "UniMod:35"

@report
def testRNaseDB():
    """
    @tests: RNaseDB
        const DigestionEnzymeRNA* getEnzyme(const String& name) nogil except +
        const DigestionEnzymeRNA* getEnzymeByRegEx(const String& cleavage_regex) nogil except +
        void getAllNames(libcpp_vector[ String ]& all_names) nogil except +
        bool hasEnzyme(const String& name) nogil except +
        bool hasRegEx(const String& cleavage_regex) nogil except +
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

    # Handling of different Unicode Strings:
    # - unicode is natively stored in OpenMS::String
    # - encoded bytesequences for utf8, utf16 and iso8859 can be stored as
    #   char arrays in OpenMS::String (and be accessed using c_str())
    # - encoded bytesequences for anything other than utf8 cannot use
    #   "toString()" as this function expects utf8
    ustr = u"bläh"
    pystr = pyopenms.String(ustr)
    assert (pystr.toString() == u"bläh")
    pystr = pyopenms.String(ustr.encode("utf8"))
    assert (pystr.toString() == u"bläh")

    pystr = pyopenms.String(ustr.encode("iso8859_15"))
    assert (pystr.c_str().decode("iso8859_15") == u"bläh")
    pystr = pyopenms.String(ustr.encode("utf16"))
    assert (pystr.c_str().decode("utf16") == u"bläh")

    # toString will throw as its not UTF8
    pystr = pyopenms.String(ustr.encode("iso8859_15"))
    didThrow = False
    try:
        pystr.toString()
    except UnicodeDecodeError:
        didThrow = True
    assert didThrow

    # toString will throw as its not UTF8
    pystr = pyopenms.String(ustr.encode("utf16"))
    didThrow = False
    try:
        pystr.toString()
    except UnicodeDecodeError:
        didThrow = True
    assert didThrow
    
    # Handling of automatic conversions of String return values
    #  -- return a native str when utf8 is used
    #  -- return a OpenMS::String object when encoding with utf8 is not possible
    ustr = u"bläh"
    s = pyopenms.MSSpectrum()
    s.setNativeID(ustr)
    r = s.getNativeID()
    # assert( isinstance(r, str) ) # native, returns str
    assert(r == u"bläh")

    s.setNativeID(ustr.encode("utf8"))
    r = s.getNativeID()
    # assert( isinstance(r, str) )
    assert(r == u"bläh")
    s.setNativeID(ustr.encode("utf16"))
    r = s.getNativeID()
    # assert( isinstance(r, bytes) )
    # assert(r.c_str().decode("utf16") == u"bläh")
    s.setNativeID(ustr.encode("iso8859_15"))
    r = s.getNativeID()
    # assert( isinstance(r, bytes) )
    assert(r.decode("iso8859_15") == u"bläh")

@report
def testGNPSExport():
    cm = pyopenms.ConsensusMap()
    
    for mz, rt, ion, linked_groups in [(222.08, 62.0, "[M+H]+", ["1","2","3"]),
                                        (244.08, 62.0, "[M+Na]+", ["1","2"]),
                                        (204.08, 62.0, "[M-H-O]+", ["3"]),
                                        (294.1, 62.0, "[M+H]+", ["4","5"])]:
        f = pyopenms.ConsensusFeature()
        f.setMZ(mz)
        f.setRT(rt)
        f.setCharge(1)
        f.setQuality(2.0)
        f.setMetaValue("best ion", ion)
        f.setMetaValue("LinkedGroups", linked_groups)
        cm.push_back(f)
    cm.setUniqueIds()

    pyopenms.IonIdentityMolecularNetworking.annotateConsensusMap(cm)

    pyopenms.IonIdentityMolecularNetworking.writeSupplementaryPairTable(cm, "SupplementaryPairsTable.csv")

    with open("SupplementaryPairsTable.csv", "r") as f:
        assert f.read() == """ID1,ID2,EdgeType,Score,Annotation
1,2,MS1 annotation,1,[M+H]+ [M+Na]+ dm/z=22.0
1,3,MS1 annotation,1,[M+H]+ [M-H-O]+ dm/z=18.0
"""
    os.remove("SupplementaryPairsTable.csv")

    pyopenms.GNPSQuantificationFile().store(cm, "FeatureQuantificationTable.txt")
    with open("FeatureQuantificationTable.txt", "r") as f:
        assert f.read() == """#MAP	id	filename	label	size
#CONSENSUS	rt_cf	mz_cf	intensity_cf	charge_cf	width_cf	quality_cf	row ID	best ion	partners	annotation network number
CONSENSUS	62.0	222.080000000000013	0.0	1	0.0	2.0	1	[M+H]+	2;3	1
CONSENSUS	62.0	244.080000000000013	0.0	1	0.0	2.0	2	[M+Na]+	1	1
CONSENSUS	62.0	204.080000000000013	0.0	1	0.0	2.0	3	[M-H-O]+	1	1
CONSENSUS	62.0	294.100000000000023	0.0	1	0.0	2.0	4	[M+H]+		2
"""
    os.remove("FeatureQuantificationTable.txt")

    # add mandatory file descriptions
    file_descriptions = {}
    for i, filename in enumerate(["1.mzML", "2.mzML"]):
        file_description = pyopenms.ColumnHeader()
        file_description.filename = filename
        file_descriptions[i] = file_description
    cm.setColumnHeaders(file_descriptions)

    pyopenms.GNPSMetaValueFile().store(cm, "MetaValueTable.tsv")
    with open("MetaValueTable.tsv", "r") as f:
        assert f.read() == """	filename	ATTRIBUTE_MAPID
0	1.mzML	MAP0
1	2.mzML	MAP1
"""
    os.remove("MetaValueTable.tsv")

