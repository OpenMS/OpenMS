#!/usr/bin/env python
# -*- coding: utf-8  -*-
from __future__ import print_function

import pyopenms
import copy
import os

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
def testElementTutorial():
    from pyopenms import *
    edb = ElementDB()

    assert edb.hasElement(b"O")
    assert edb.hasElement(b"S")

    oxygen = edb.getElement(b"O")
    assert oxygen.getName() == b'Oxygen'
    assert oxygen.getSymbol() == b"O"
    assert abs(oxygen.getMonoWeight() - 15.994915) < 1e-5
    isotopes = oxygen.getIsotopeDistribution()

    sulfur = edb.getElement("S")
    sulfur.getName() 
    sulfur.getSymbol()
    sulfur.getMonoWeight()
    isotopes = sulfur.getIsotopeDistribution()
    exp = [(31.97207073, 0.9493), (32.971458, 0.0076), (33.967867, 0.0429), (35.967081, 0.0002)]
    for e, iso in zip(exp, isotopes.getContainer()):
        assert abs(e[0] - iso.getMZ()) < 1e-5
        assert abs(e[1] - iso.getIntensity()) < 1e-5

@report
def testEmpiricalFormulaTutorial():

    from pyopenms import *

    methanol = EmpiricalFormula("CH3OH")
    water = EmpiricalFormula("H2O")
    wm = EmpiricalFormula(str(water) + str(methanol))
    print(wm)

    isotopes = wm.getIsotopeDistribution( CoarseIsotopePatternGenerator(3) )
    for iso in isotopes.getContainer():
        print (iso)

    # .. Wait for pyOpenMS 2.4
    wm = water + methanol # only in pyOpenMS 2.4
    m = wm.getElementalComposition()
    assert m["C"] == 1
    assert m["H"] == 6
    assert m["O"] == 2

@report
def testResidueTutorial():

    from pyopenms import *
    lys = ResidueDB().getResidue("Lysine")
    assert lys.getName() == b"Lysine"
    lys.getThreeLetterCode() == b"LYS"
    lys.getOneLetterCode() == b"K"
    lys.getAverageWeight(Residue.ResidueType.Full)
    assert abs(lys.getMonoWeight(Residue.ResidueType.Full) - 146.1055284466) < 1e-5
    lys.getPka()

@report
def testAASequenceTutorial():

    from pyopenms import *
    seq = AASequence.fromString("DFPIANGER")
    prefix = seq.getPrefix(4)
    suffix = seq.getSuffix(5)
    concat = seq + seq

    print(seq)
    print(concat)
    print(suffix)
    seq.getMonoWeight(Residue.ResidueType.Full, 0)
    seq.getMonoWeight(Residue.ResidueType.Full, 2) / 2.0
    assert abs(concat.getMonoWeight(Residue.ResidueType.Full, 0)- 2016.9653632108004) < 1e-5

    seq_formula = seq.getFormula(Residue.ResidueType.Full, 0)
    print(seq_formula)

    isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(6) )
    for iso in isotopes.getContainer():
        print (iso)

    assert abs(isotopes.getContainer()[0].getIntensity() - 0.5681651355598922) < 1e-5
    suffix = seq.getSuffix(3) # y3 ion
    y3_formula = suffix.getFormula(Residue.ResidueType.YIon, 2) # y3++ ion
    suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0
    suffix.getMonoWeight(Residue.ResidueType.XIon, 2) / 2.0 # ATTENTION
    suffix.getMonoWeight(Residue.ResidueType.BIon, 2) / 2.0 # ATTENTION
    assert str(y3_formula) == b"C13H24N6O6", str(y3_formula)
    assert str(seq_formula) == b"C44H67N13O15"


    from pyopenms import *
    seq = AASequence.fromString("PEPTIDESEKUEM(Oxidation)CER")
    print(seq.toString())
    print(seq.toUnmodifiedString())
    print(seq.toBracketString(True, []))
    print(seq.toBracketString(False, []))
    print(seq.toUniModString()) # with 2.4

    assert seq.toBracketString(True, []) == b"PEPTIDESEKUEM[147]CER"

@report
def testTheoreticalSpectrumGenTutorial():

    from pyopenms import *

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

    from pyopenms import *
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

    from pyopenms import *
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

    from pyopenms import *
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

    from pyopenms import *

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


