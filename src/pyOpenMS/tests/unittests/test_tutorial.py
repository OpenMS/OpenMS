#!/usr/bin/env python3
# -*- coding: utf-8  -*-
import copy
import os

from pyopenms import String as s
from pyopenms import *

try:
    long
except NameError:
    long = int

from functools import wraps


def report(f):
    @wraps(f)
    def wrapper(*a, **kw):
        print("run b", f.__name__)
        f(*a, **kw)
    return wrapper

@report
def testElementTutorial():

    edb = ElementDB()

    assert edb.hasElement("O")
    assert edb.hasElement("S")

    oxygen = edb.getElement("O")
    assert oxygen.getName() == "Oxygen"
    assert oxygen.getSymbol() == "O"
    oxygen.getMonoWeight()
    isotopes = oxygen.getIsotopeDistribution()

    sulfur = edb.getElement("S")
    sulfur.getName()
    sulfur.getSymbol()
    sulfur.getMonoWeight()
    isotopes = sulfur.getIsotopeDistribution()
    for iso in isotopes.getContainer():
        print (iso.getMZ(), ":", iso.getIntensity())

@report
def testEmpiricalFormulaTutorial():

    print("testEmpiricalFormulaTutorial")
    methanol = EmpiricalFormula("CH3OH")
    water = EmpiricalFormula("H2O")
    ethanol = EmpiricalFormula("CH2" + methanol.toString())
    wm = water + methanol
    print(wm)
    print(wm.getElementalComposition())

    wm = water + methanol # only in pyOpenMS 2.4
    m = wm.getElementalComposition()
    assert m["C"] == 1
    assert m["H"] == 6
    assert m["O"] == 2

    wm = EmpiricalFormula("CH3OH") + EmpiricalFormula("H2O")

    isotopes = wm.getIsotopeDistribution( CoarseIsotopePatternGenerator(3) )
    for iso in isotopes.getContainer():
        print (iso.getMZ(), ":", iso.getIntensity())

    isotopes = wm.getIsotopeDistribution( CoarseIsotopePatternGenerator(3, True) )
    for iso in isotopes.getContainer():
        print (iso.getMZ(), ":", iso.getIntensity())

@report
def testResidueTutorial():

    lys = ResidueDB().getResidue("Lysine")
    assert lys.getName() == "Lysine"
    lys.getThreeLetterCode() == "LYS"
    lys.getOneLetterCode() == "K"
    lys.getAverageWeight(Residue.ResidueType.Full)
    assert abs(lys.getMonoWeight(Residue.ResidueType.Full) - 146.1055284466) < 1e-5
    lys.getPka()

@report
def testAASequenceLen():
    """Test the __len__() method implementation for AASequence"""
    
    # Test basic len() functionality
    seq = AASequence.fromString("PEPTIDE")
    assert len(seq) == 7
    assert len(seq) == seq.size()
    
    # Test with different sequences
    short_seq = AASequence.fromString("PEP")
    assert len(short_seq) == 3
    
    long_seq = AASequence.fromString("PEPTIDESEQUENCE")
    assert len(long_seq) == 15
    
    # Test with modified sequences
    modified_seq = AASequence.fromString("PEPTIDEM(Oxidation)")
    assert len(modified_seq) == 8  # Modification doesn't change sequence length
    
    # Test with empty sequence (if possible)
    empty_seq = AASequence.fromString("")
    assert len(empty_seq) == 0
    
    print("All len() tests passed!")

@report
def testMSExperimentLen():
    """Test the __len__() method implementation for MSExperiment"""
    
    # Test with empty experiment
    exp = MSExperiment()
    assert len(exp) == 0
    assert len(exp) == exp.size()
    
    # Test with spectra added
    spectrum1 = MSSpectrum()
    spectrum1.setRT(1.0)
    exp.addSpectrum(spectrum1)
    assert len(exp) == 1
    assert len(exp) == exp.size()
    
    spectrum2 = MSSpectrum()
    spectrum2.setRT(2.0)
    exp.addSpectrum(spectrum2)
    assert len(exp) == 2
    assert len(exp) == exp.size()
    
    # Test with multiple spectra
    exp2 = MSExperiment()
    for i in range(5):
        spectrum = MSSpectrum()
        spectrum.setRT(float(i))
        exp2.addSpectrum(spectrum)
    assert len(exp2) == 5
    assert len(exp2) == exp2.size()
    assert len(exp2) == exp2.getNrSpectra()
    
    print("All MSExperiment len() tests passed!")

@report
def testAASequenceTutorial():

    seq = AASequence.fromString("DFPIANGER")
    # Test the new len() functionality
    assert len(seq) == 9  # "DFPIANGER" has 9 amino acids
    assert len(seq) == seq.size()  # len() should return same as size()
    
    prefix = seq.getPrefix(4)
    suffix = seq.getSuffix(5)
    concat = seq + seq
    
    # Test len() with different sequences
    assert len(prefix) == 4
    assert len(suffix) == 5
    assert len(concat) == 18  # concatenated sequence should be 2 * 9

    print(seq)
    print(concat)
    print(suffix)
    print("len(seq) =", len(seq), ", seq.size() =", seq.size())
    seq.getMonoWeight() # weight of M
    seq.getMonoWeight(Residue.ResidueType.Full, 2) # weight of M+2H
    mz = seq.getMonoWeight(Residue.ResidueType.Full, 2) / 2.0 # m/z of M+2H
    concat.getMonoWeight()

    print("Monoisotopic m/z of (M+2H)2+ is", mz)

    seq_formula = seq.getFormula()
    print("Peptide", seq, "has molecular formula", seq_formula)
    print("="*35)

    isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(6) )
    for iso in isotopes.getContainer():
            print ("Isotope", iso.getMZ(), ":", iso.getIntensity())

    suffix = seq.getSuffix(3) # y3 ion "GER"
    print("="*35)
    print("y3 ion :", suffix)
    y3_formula = suffix.getFormula(Residue.ResidueType.YIon, 2) # y3++ ion
    suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 # CORRECT
    suffix.getMonoWeight(Residue.ResidueType.XIon, 2) / 2.0 # CORRECT
    suffix.getMonoWeight(Residue.ResidueType.BIon, 2) / 2.0 # INCORRECT
    assert str(y3_formula) == "C13H24N6O6", str(y3_formula)
    assert str(seq_formula) == "C44H67N13O15"


    print("y3 mz :", suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 )
    print(y3_formula)
    print(seq_formula)

    seq = AASequence.fromString("PEPTIDESEKUEM(Oxidation)CER")
    print(seq)
    print(seq.toString())
    print(seq.toUnmodifiedString())
    print(seq.toUniModString())
    print(seq.toBracketString())
    print(seq.toBracketString(False))

    print(AASequence.fromString("DFPIAM(UniMod:35)GER"))
    print(AASequence.fromString("DFPIAM[+16]GER"))
    print(AASequence.fromString("DFPIAM[+15.99]GER"))
    print(AASequence.fromString("DFPIAM[147]GER"))
    print(AASequence.fromString("DFPIAM[147.035405]GER"))

    s = AASequence.fromString(".(Dimethyl)DFPIAMGER.")
    print(s, s.hasNTerminalModification())
    s = AASequence.fromString(".DFPIAMGER.(Label:18O(2))")
    print(s, s.hasCTerminalModification())
    s = AASequence.fromString(".DFPIAMGER(Phospho).")
    print(s, s.hasCTerminalModification())


    bsa = FASTAEntry()
    bsa.sequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
    bsa.description = "BSA Bovine Albumin (partial sequence)"
    bsa.identifier = "BSA"
    alb = FASTAEntry()
    alb.sequence = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGE"
    alb.description = "ALB Human Albumin (partial sequence)"
    alb.identifier = "ALB"

    entries = [bsa, alb]

    f = FASTAFile()
    f.store("example.fasta", entries)

    entries = []
    f = FASTAFile()
    f.load("example.fasta", entries)
    print( len(entries) )
    for e in entries:
      print (e.identifier, e.sequence)

@report
def testTheoreticalSpectrumGenTutorial():


    tsg = TheoreticalSpectrumGenerator()
    spec1 = MSSpectrum()
    spec2 = MSSpectrum()
    peptide = AASequence.fromString("DFPIANGER")
    # standard behavior is adding b- and y-ions of charge 1
    p = Param()
    p.setValue("add_b_ions", "false", "Add peaks of b-ions to the spectrum")
    tsg.setParameters(p)
    tsg.getSpectrum(spec1, peptide, 1, 1)
    p.setValue("add_b_ions", "true", "Add peaks of a-ions to the spectrum")
    p.setValue("add_metainfo", "true", "")
    tsg.setParameters(p)
    tsg.getSpectrum(spec2, peptide, 1, 2)
    print("Spectrum 1 has", spec1.size(), "peaks.")
    print("Spectrum 2 has", spec2.size(), "peaks.")

    # Iterate over annotated ions and their masses
    for ion, peak in zip(spec2.getStringDataArrays()[0], spec2):
        print(ion, peak.getMZ())



@report
def testDigestionTutorial():

    print("testDigestionTutorial")
    dig = ProteaseDigestion()
    dig.getEnzymeName() # Trypsin
    bsa = b'DEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA'
    bsa = AASequence.fromString(bsa)
    result = []
    dig.digest(bsa, result)
    print(result[4])
    len(result) # 74 peptides
    assert len(result) == 74

    names = []
    ProteaseDB().getAllNames(names)
    len(names) # 25 by default
    e = ProteaseDB().getEnzyme('Lys-C')
    e.getRegExDescription()
    e.getRegEx()

    dig = ProteaseDigestion()
    dig.setEnzyme('Lys-C')
    result = []
    dig.digest(bsa, result)
    print(result[4])
    len(result) # 53 peptides
    assert len(result) == 53

def gaussian(x, mu, sig):
    import numpy as np
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

@report
def testDatastructuresTutorial():

    ###################################
    # Spectrum
    # ********
    ###################################

    spectrum = MSSpectrum()
    mz = range(1500, 500, -100)
    i = [0 for mass in mz]
    spectrum.set_peaks([mz, i])

    # Sort the peaks according to ascending mass-to-charge ratio
    spectrum.sortByPosition()

    # Iterate over spectrum of those peaks
    for p in spectrum:
        print(p.getMZ(), p.getIntensity())

    # More efficient peak access
    for mz, i in zip(*spectrum.get_peaks()):
        print(mz, i)

    # Access a peak by index
    print(spectrum[1].getMZ(), spectrum[1].getIntensity())


    ###################################
    ###################################

    spectrum = MSSpectrum()
    spectrum.setDriftTime(25) # 25 ms
    spectrum.setRT(205.2) # 205.2 s
    spectrum.setMSLevel(3) # MS3
    p = Precursor()
    p.setIsolationWindowLowerOffset(1.5)
    p.setIsolationWindowUpperOffset(1.5) 
    p.setMZ(600) # isolation at 600 +/- 1.5 Th
    p.setActivationEnergy(40) # 40 eV
    p.setCharge(4) # 4+ ion
    spectrum.setPrecursors( [p] )

    fda = FloatDataArray()
    fda.setName("Signal to Noise Array")
    fda.push_back(15)
    sda = StringDataArray()
    sda.setName("Peak annotation")
    sda.push_back("y15++")
    spectrum.setFloatDataArrays( [fda] )
    spectrum.setStringDataArrays( [sda] )
    spectrum.set_peaks( ([401.5], [900]) )

    # Store as mzML
    exp = MSExperiment()
    exp.addSpectrum(spectrum)
    spectrum = MSSpectrum()
    spectrum.set_peaks( ([1, 2], [1, 2]) )
    exp.addSpectrum(spectrum)
    MzMLFile().store("testfile.mzML", exp)

    ###################################
    # Peak Map
    # *********
    ###################################

    # The following examples creates a MSExperiment containing four MSSpectrum instances.
    exp = MSExperiment()
    for i in range(6):
        spectrum = MSSpectrum()
        spectrum.setRT(i)
        spectrum.setMSLevel(1)
        for mz in range(500, 900, 100):
          peak = Peak1D()
          peak.setMZ(mz + i)
          peak.setIntensity(100 - 25*abs(i-2.5) )
          spectrum.push_back(peak)
        exp.addSpectrum(spectrum)

    # Iterate over spectra
    for s in exp:
        for p in s:
            print (s.getRT(), p.getMZ(), p.getIntensity())

    # Sum intensity of all spectra between RT 2.0 and 3.0
    print(sum([p.getIntensity() for s in exp if s.getRT() >= 2.0 and s.getRT() <= 3.0 for p in s]))

    # Store as mzML
    MzMLFile().store("testfile2.mzML", exp)

    ###################################
    # Chromatogram
    # ************
    ###################################

    chromatogram = MSChromatogram()
    rt = range(1500, 500, -100)
    i = [gaussian(rtime, 1000, 150) for rtime in rt]
    chromatogram.set_peaks([rt, i])

    # Sort the peaks according to ascending retention time
    chromatogram.sortByPosition()

    # Iterate over chromatogram of those peaks
    for p in chromatogram:
        print(p.getRT(), p.getIntensity())

    # Access a peak by index
    print(chromatogram[1].getRT(), chromatogram[1].getIntensity())

    chromatogram.setNativeID("Trace XIC@405.2")
    p = Precursor()
    p.setIsolationWindowLowerOffset(1.5)
    p.setIsolationWindowUpperOffset(1.5) 
    p.setMZ(405.2) # isolation at 405.2 +/- 1.5 Th
    p.setActivationEnergy(40) # 40 eV
    p.setCharge(2) # 2+ ion
    p.setMetaValue("description", "test_description")
    p.setMetaValue("description", chromatogram.getNativeID())
    p.setMetaValue("peptide_sequence", chromatogram.getNativeID())
    chromatogram.setPrecursor(p)
    p = Product()
    p.setMZ(603.4) # transition from 405.2 -> 603.4
    chromatogram.setProduct(p)

    # Store as mzML
    exp = MSExperiment()
    exp.addChromatogram(chromatogram)
    MzMLFile().store("testfile3.mzML", exp)

@report
def testMSSpectrumLen():
    """Test the __len__() method implementation for MSSpectrum"""
    
    # Test with empty spectrum
    spectrum = MSSpectrum()
    assert len(spectrum) == 0
    assert len(spectrum) == spectrum.size()
    
    # Test with peaks added
    p1 = Peak1D()
    p1.setMZ(500.0)
    p1.setIntensity(1e5)
    spectrum.push_back(p1)
    assert len(spectrum) == 1
    assert len(spectrum) == spectrum.size()
    
    p2 = Peak1D()
    p2.setMZ(600.0)
    p2.setIntensity(2e5)
    spectrum.push_back(p2)
    assert len(spectrum) == 2
    assert len(spectrum) == spectrum.size()
    
    # Test with multiple peaks using set_peaks
    spectrum2 = MSSpectrum()
    mz_values = [100.0, 200.0, 300.0, 400.0, 500.0]
    intensity_values = [1.0, 2.0, 3.0, 4.0, 5.0]
    spectrum2.set_peaks([mz_values, intensity_values])
    assert len(spectrum2) == 5
    assert len(spectrum2) == spectrum2.size()
    
    print("All MSSpectrum len() tests passed!")

@report
def testMSChromatogramLen():
    """Test the __len__() method implementation for MSChromatogram"""
    
    # Test with empty chromatogram
    chromatogram = MSChromatogram()
    assert len(chromatogram) == 0
    assert len(chromatogram) == chromatogram.size()
    
    # Test with peaks added
    p1 = ChromatogramPeak()
    p1.setRT(1.0)
    p1.setIntensity(100.0)
    chromatogram.push_back(p1)
    assert len(chromatogram) == 1
    assert len(chromatogram) == chromatogram.size()
    
    p2 = ChromatogramPeak()
    p2.setRT(2.0)
    p2.setIntensity(200.0)
    chromatogram.push_back(p2)
    assert len(chromatogram) == 2
    assert len(chromatogram) == chromatogram.size()
    
    # Test with multiple peaks using set_peaks
    chromatogram2 = MSChromatogram()
    rt_values = [1.0, 2.0, 3.0, 4.0, 5.0]
    intensity_values = [10.0, 20.0, 30.0, 40.0, 50.0]
    chromatogram2.set_peaks([rt_values, intensity_values])
    assert len(chromatogram2) == 5
    assert len(chromatogram2) == chromatogram2.size()
    
    print("All MSChromatogram len() tests passed!")
