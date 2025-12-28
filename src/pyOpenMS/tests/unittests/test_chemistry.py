#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Chemistry-related tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms
import copy
import numpy as np

# Import shared test helper functions
from test_helpers import (
  report, 
  _testMetaInfoInterface,
  _testUniqueIdInterface,
  _testStrOutput
)

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
   AASequence.getAAFrequencies
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

  # constructor from C string using the static method
  seq = pyopenms.AASequence.fromString("PEPTIDESEKUEM(Oxidation)CER")
  assert seq.toString() == "PEPTIDESEKUEM(Oxidation)CER"
  assert seq.toUnmodifiedString() == "PEPTIDESEKUEMCER"
  assert seq.toBracketString() == "PEPTIDESEKUEM[147]CER"
  assert seq.toBracketString(True) == "PEPTIDESEKUEM[147]CER"

  # constructor from String
  seq2 = pyopenms.AASequence("PEPTIDESEKUEM(Oxidation)CER")
  assert seq2.toString() == "PEPTIDESEKUEM(Oxidation)CER"
  assert seq2.toUnmodifiedString() == "PEPTIDESEKUEMCER"
  assert seq == seq2
  assert seq2.toBracketString() == "PEPTIDESEKUEM[147]CER"
  assert seq2.toBracketString(True) == "PEPTIDESEKUEM[147]CER"

  # constructor from String (Permissive)
  seq3 = pyopenms.AASequence("PEPTIDE#SEKUEM(Oxidation)CER", True)
  assert seq3.toString() == "PEPTIDEXSEKUEM(Oxidation)CER"
  assert seq3.toUnmodifiedString() == "PEPTIDEXSEKUEMCER"
  assert seq3.toBracketString() == "PEPTIDEXSEKUEM[147]CER"
  assert seq3.toBracketString(True) == "PEPTIDEXSEKUEM[147]CER"

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

  # Test __repr__ and __str__ methods
  aas_repr = pyopenms.AASequence.fromString("PEPTM(Oxidation)IDE")
  repr_str = repr(aas_repr)
  assert "AASequence(" in repr_str
  assert "sequence=" in repr_str
  assert "PEPTM(Oxidation)IDE" in repr_str
  assert "length=" in repr_str
  assert "mono_mass=" in repr_str
  assert "modified=True" in repr_str
  str_str = str(aas_repr)
  assert str_str == repr_str

  # Test unmodified sequence
  aas_unmod = pyopenms.AASequence.fromString("PEPTIDE")
  repr_unmod = repr(aas_unmod)
  assert "modified=True" not in repr_unmod

  # Test getAAFrequencies - matching C++ test case
  aas_freq = pyopenms.AASequence.fromString("THREEAAAWITHYYY")
  freq_table = {}
  aas_freq.getAAFrequencies(freq_table)
  assert freq_table[b"T"] == 2
  assert freq_table[b"H"] == 2
  assert freq_table[b"R"] == 1
  assert freq_table[b"E"] == 2
  assert freq_table[b"A"] == 3
  assert freq_table[b"W"] == 1
  assert freq_table[b"I"] == 1
  assert freq_table[b"Y"] == 3
  assert len(freq_table) == 8


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
  oms_string = pyopenms.String("blu")
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
  assert pyopenms.String(e.getSymbol()) == pyopenms.String(u"blü")
  assert pyopenms.String(e.getSymbol()).toString() == u"blü"

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
def testResidueRepr():
  """
  @tests: Residue.__repr__
  """
  # Get a residue from the database
  rdb = pyopenms.ResidueDB()
  glycine = rdb.getResidue(pyopenms.String("Glycine"))

  # Test __repr__ method
  repr_str = repr(glycine)
  assert "Residue(" in repr_str
  assert "name=" in repr_str
  assert "Glycine" in repr_str
  assert "one_letter=" in repr_str
  assert "'G'" in repr_str
  assert "three_letter=" in repr_str
  assert "'Gly'" in repr_str
  assert "formula=" in repr_str
  assert "mono_mass=" in repr_str

  # Test __str__ method for unmodified residue
  str_str = str(glycine)
  assert str_str == "G"

  # Test with modified residue - get oxidized methionine
  methionine = rdb.getResidue(pyopenms.String("Methionine"))
  str_str = str(methionine)
  assert str_str == "M"
  repr_str = repr(methionine)
  assert "Residue(" in repr_str
  assert "'M'" in repr_str
  assert "'Met'" in repr_str


@report
def testResidueModificationRepr():
  """
  @tests: ResidueModification.__repr__
  """
  # Get a modification from the database
  mdb = pyopenms.ModificationsDB()
  mod = mdb.getModification("Oxidation", "M", pyopenms.ResidueModification.TermSpecificity.ANYWHERE)

  # Test __repr__ and __str__ methods
  repr_str = repr(mod)
  assert "ResidueModification(" in repr_str
  str_str = str(mod)
  assert str_str == repr_str


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

  # Test __repr__ and __str__ methods
  iso_repr = pyopenms.IsotopeDistribution()
  iso_repr.insert(100.0, 0.9)
  iso_repr.insert(101.0, 0.08)
  iso_repr.insert(102.0, 0.02)

  repr_str = repr(iso_repr)
  assert "IsotopeDistribution(" in repr_str
  assert "num_isotopes=" in repr_str
  assert "mass_range=" in repr_str
  str_str = str(iso_repr)
  assert str_str == repr_str


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

  # Test __repr__ and __str__ methods
  ef_repr = pyopenms.EmpiricalFormula("C6H12O6")
  repr_str = repr(ef_repr)
  assert "EmpiricalFormula(" in repr_str
  assert "formula=" in repr_str
  assert "C6H12O6" in repr_str
  assert "mono_mass=" in repr_str
  str_str = str(ef_repr)
  # __str__ returns just the formula, __repr__ returns detailed info
  assert str_str == "C6H12O6"
  assert str_str != repr_str


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
def testAdduct():
  """
  @tests: Adduct
   Adduct.__init__
  """
  a = pyopenms.Adduct()

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
  # const Map[ String, Element * ]  getNames() except + nogil 
  # const Map[ String, Element * ] getSymbols() except + nogil 
  # const Map[unsigned int, Element * ] getAtomicNumbers() except + nogil 


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

  assert rdb.hasResidue(pyopenms.String("Glycine"))
  glycine = rdb.getResidue(pyopenms.String("Glycine"))

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
  mdb.searchModifications(mods, pyopenms.String("Phosphorylation"), pyopenms.String("T"), pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
  assert len(mods) == 1

  mods = set([])
  mdb.searchModifications(mods, pyopenms.String("NIC"), pyopenms.String("T"), pyopenms.ResidueModification.TermSpecificity.N_TERM)
  assert len(mods) == 1

  mods = set([])
  mdb.searchModifications(mods, pyopenms.String("NIC"), pyopenms.String("T"), pyopenms.ResidueModification.TermSpecificity.N_TERM)
  assert len(mods) == 1

  mods = set([])
  mdb.searchModifications(mods, pyopenms.String("Acetyl"), pyopenms.String("T"), pyopenms.ResidueModification.TermSpecificity.N_TERM)
  assert len(mods) == 1
  assert list(mods)[0].getFullId() == "Acetyl (N-term)"

  m = mdb.getModification(pyopenms.String("Carboxymethyl (C)"), "", pyopenms.ResidueModification.TermSpecificity.NUMBER_OF_TERM_SPECIFICITY)
  assert m.getFullId() == "Carboxymethyl (C)"

  m = mdb.getModification( pyopenms.String("Phosphorylation"), pyopenms.String("S"), pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
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
def testProteaseDB():
  edb = pyopenms.ProteaseDB()

  f = pyopenms.EmpiricalFormula()
  synonyms = set(["dummy", "other"])

  assert edb.hasEnzyme(pyopenms.String("Trypsin"))

  trypsin = edb.getEnzyme(pyopenms.String("Trypsin"))

  names = []
  edb.getAllNames(names)
  assert b"Trypsin" in names


if __name__ == "__main__":
  # Run all tests in this module
  import sys
  module = sys.modules[__name__]
    
  test_functions = [getattr(module, name) for name in dir(module) 
           if name.startswith('test') and callable(getattr(module, name))]
    
  print(f"Running {len(test_functions)} chemistry tests...")
  for test_func in test_functions:
    try:
      test_func()
      print(f"✅ {test_func.__name__}")
    except Exception as e:
      print(f"❌ {test_func.__name__}: {e}")
