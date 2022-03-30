#!/usr/bin/env python
# -*- coding: utf-8  -*-
from __future__ import print_function

import pyopenms
import copy
import os

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

    # replace oxygen
    e2 = edb.addElement(b"Oxygen", b"O", 8, {16 : 0.7, 19 : 0.3}, {16 : 16.01, 19 : 19.01}, True)
    e2 = edb.getElement(pyopenms.String("O"))
    assert e2.getName() == "Oxygen"
    assert e2.getIsotopeDistribution()
    assert len(e2.getIsotopeDistribution().getContainer()) == 2
    assert abs(e2.getIsotopeDistribution().getContainer()[1].getIntensity() - 0.3) < 1e-5

    iso = pyopenms.ElementDB().getIsotope("(13)C")
    assert iso.getNeutrons() == 7
    assert iso.isStable()

    iso = pyopenms.ElementDB().getIsotope("(238)U")
    assert iso.getNeutrons() == 146
    assert not iso.isStable()
    assert iso.getSymbol
    assert iso.getSymbol() == "(238)U"
    assert abs(iso.getHalfLife() - 1.40996461536e+17) < 1e-6

    ele = pyopenms.ElementDB().getElement("U")
    iso = [i for i in ele.getIsotopes() if i.getNeutrons() ==  146]
    assert len(iso) == 1
    assert iso[0].getElement().getSymbol() == "U"

    edb.addIsotope(b"NewElement", b"NE", 300, 0.5, 404, 100, 0, False)

    #  not yet implemented
    #
    # const Map[ String, Element * ]  getNames() nogil except +
    # const Map[ String, Element * ] getSymbols() nogil except +
    # const Map[unsigned int, Element * ] getAtomicNumbers() nogil except +


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
    assert e.getRegEx() == u'(?<=G)'
    assert e.getThreePrimeGain() == u'p'

    assert db.hasRegEx(u'(?<=G)')
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
    # assert abs(seq.getMonoWeight(pyopenms.Residue.ResidueType.Full, 0) - 1952.7200317517998) < 1e-5
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

