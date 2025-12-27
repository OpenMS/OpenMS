#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Protein/Peptide Identification tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms

# Import shared test helper functions
from test_helpers import report, _testMetaInfoInterface, _testParam, _testStrOutput, _testUniqueIdInterface


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

    # Test __repr__ and __str__ methods
    ph_repr = pyopenms.PeptideHit()
    ph_repr.setScore(18.1)
    ph_repr.setRank(1)
    ph_repr.setCharge(2)
    ph_repr.setSequence(pyopenms.AASequence.fromString("PEPTIDER"))
    
    # Basic repr test
    repr_str = repr(ph_repr)
    assert "PeptideHit(" in repr_str
    assert "score=" in repr_str
    assert "sequence=" in repr_str
    assert "charge=" in repr_str
    
    # Test with protein evidences
    pe1 = pyopenms.PeptideEvidence()
    pe1.setProteinAccession('PH_6057')
    pe1.setStart(71)
    pe1.setEnd(80)
    pe1.setAABefore(b'R')
    pe1.setAAAfter(b'N')
    
    ph_repr.addPeptideEvidence(pe1)
    repr_str = repr(ph_repr)
    assert "evidences=" in repr_str
    assert "PH_6057" in repr_str
    assert "PeptideEvidence(" in repr_str
    
    # Test str method
    str_str = str(ph_repr)
    assert str_str == repr_str


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

    # Test __repr__ and __str__ methods
    pe_repr = pyopenms.PeptideEvidence()
    pe_repr.setProteinAccession('PH_6057')
    pe_repr.setStart(71)
    pe_repr.setEnd(80)
    pe_repr.setAABefore(b'R')
    pe_repr.setAAAfter(b'N')
    
    repr_str = repr(pe_repr)
    assert "PeptideEvidence(" in repr_str
    assert "protein=" in repr_str
    assert "PH_6057" in repr_str
    assert "start=" in repr_str
    assert "end=" in repr_str
    assert "aa_before=" in repr_str
    assert "aa_after=" in repr_str
    
    # Test str method
    str_str = str(pe_repr)
    assert str_str == repr_str


@report
def testPeptideIdentification():
    """
    @tests: PeptideIdentification
     PeptideIdentification.__init__
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

    assert isinstance(pi.getSignificanceThreshold(), float)
    _testStrOutput(pi.getScoreType())
    pi.setScoreType("A")
    assert isinstance(pi.isHigherScoreBetter(), int)
    _testStrOutput(pi.getIdentifier())
    pi.setIdentifier("id")
    pi.sort()
    assert not pi.empty()

    pi.setSignificanceThreshold(6.0)

    # Test __repr__ and __str__ methods
    pi_repr = pyopenms.PeptideIdentification()
    pi_repr.setRT(1234.5)
    pi_repr.setMZ(445.678)
    pi_repr.setScoreType("XTandem")
    
    hit = pyopenms.PeptideHit()
    hit.setScore(50.5)
    hit.setRank(1)
    hit.setSequence(pyopenms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    pi_repr.insertHit(hit)
    
    repr_str = repr(pi_repr)
    assert "PeptideIdentification(" in repr_str
    assert "rt=" in repr_str
    assert "mz=" in repr_str
    assert "score_type=" in repr_str
    assert "num_hits=" in repr_str
    assert "top_hit=" in repr_str
    assert "PEPTIDE" in repr_str
    
    # Test str method
    str_str = str(pi_repr)
    assert str_str == repr_str


@report
def testPeptideIdentificationList():
    """
    @tests: PeptideIdentificationList
     PeptideIdentificationList.__init__
     PeptideIdentificationList.__getitem__
     PeptideIdentificationList.__iter__
     PeptideIdentificationList.__len__
     PeptideIdentificationList.clear
     PeptideIdentificationList.empty
     PeptideIdentificationList.push_back
     PeptideIdentificationList.size
    """
    # Test default constructor
    pil = pyopenms.PeptideIdentificationList()
    assert pil.empty()
    assert pil.size() == 0

    # Create some PeptideIdentification objects for testing
    pi1 = pyopenms.PeptideIdentification()
    pi1.setRT(100.0)
    pi1.setMZ(200.0)
    pi1.setIdentifier("test1")
    
    pi2 = pyopenms.PeptideIdentification()
    pi2.setRT(150.0)
    pi2.setMZ(250.0)
    pi2.setIdentifier("test2")
    
    # Test push_back
    pil.push_back(pi1)
    assert pil.size() == 1
    assert not pil.empty()
    
    # Test push_back again (append is same as push_back)
    pil.push_back(pi2)
    assert pil.size() == 2
    
    # Test __len__
    assert len(pil) == 2
    
    # Test __getitem__
    retrieved_pi1 = pil[0]
    assert retrieved_pi1.getRT() == 100.0
    assert retrieved_pi1.getIdentifier() == "test1"
    
    retrieved_pi2 = pil[1]
    assert retrieved_pi2.getRT() == 150.0
    
    # Test __iter__
    rts = [pi.getRT() for pi in pil]
    assert rts == [100.0, 150.0]
    
    # Test push_back with additional element
    pi3 = pyopenms.PeptideIdentification()
    pi3.setRT(200.0)
    pil.push_back(pi3)
    assert pil.size() == 3
    
    # Test clear
    pil.clear()
    assert pil.empty()
    assert pil.size() == 0


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

    # Test __repr__ and __str__ methods
    ph_repr = pyopenms.ProteinHit()
    ph_repr.setAccession("P12345")
    ph_repr.setScore(150.5)
    ph_repr.setRank(1)
    ph_repr.setCoverage(45.2)
    ph_repr.setDescription("Example protein")
    
    repr_str = repr(ph_repr)
    assert "ProteinHit(" in repr_str
    assert "accession=" in repr_str
    assert "P12345" in repr_str
    assert "score=" in repr_str
    assert "coverage=" in repr_str
    
    # Test str method
    str_str = str(ph_repr)
    assert str_str == repr_str


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

    # Test getAllNamesOf method
    peak_mass_names = pyopenms.ProteinIdentification.getAllNamesOfPeakMassType()
    assert len(peak_mass_names) == pyopenms.ProteinIdentification.PeakMassType.SIZE_OF_PEAKMASSTYPE
    assert peak_mass_names[pyopenms.ProteinIdentification.PeakMassType.MONOISOTOPIC].decode() == "Monoisotopic"
    assert peak_mass_names[pyopenms.ProteinIdentification.PeakMassType.AVERAGE].decode() == "Average"


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
def testIDRipper():
    """
    @tests: IDRipper
     IDRipper.__init__
     IDRipper.rip
    """
    ff = pyopenms.IDRipper()

    assert pyopenms.IDRipper().rip is not None


@report
def testIDFilter():
    """
    @tests: IDFilter
     IDFilter.__init__
    """
    ff = pyopenms.IDFilter()


@report
def testIDDecoyProbability():
    """
    @tests: IDDecoyProbability
      IDDecoyProbability.__init__
    """
    ff = pyopenms.IDDecoyProbability()

    assert pyopenms.IDDecoyProbability().apply is not None


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

    target = [float(i) for i in range(10)]
    model.getGumbelGnuplotFormula(GaussFitResult) 
    model.getGaussGnuplotFormula(GaussFitResult) 
    model.getBothGnuplotFormula(GaussFitResult, GaussFitResult) 
    model.plotTargetDecoyEstimation(target, target)
    model.getSmallestScore()


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
def testPeptideIndexing():
    """
    @tests: PeptideIndexing
     PeptideIndexing.__init__
    """
    e = pyopenms.PeptideIndexing()
    assert e


@report
def testPeptideProteinResolution():
    """
    @tests: PeptideProteinResolution
     PeptideProteinResolution.__init__
    """
    e = pyopenms.PeptideProteinResolution(False)
    assert e


@report
def test_peptide_identifications_to_df():
    """
    @tests: peptide_identifications_to_df
     peptide_identifications_to_df
     update_scores_from_df
    """
    # convert to dataframe
    peps = pyopenms.PeptideIdentificationList()

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

    peps.push_back(p)

    p1 = pyopenms.PeptideIdentification()
    p1.setRT(1243.56)
    p1.setMZ(240.0)
    p1.setScoreType("ScoreType")
    p1.setHigherScoreBetter(False)
    p1.setIdentifier("IdentificationRun2")

    peps.push_back(p1)

    assert pyopenms.peptide_identifications_to_df(peps).shape == (2,12)
    assert pyopenms.peptide_identifications_to_df(peps, decode_ontology=False).shape == (2,12)
    assert pyopenms.peptide_identifications_to_df(peps)['protein_accession'][0] == 'sp|Accession1,sp|Accession2'
    assert pyopenms.peptide_identifications_to_df(peps, export_unidentified=False).shape == (1,12)

    # update from dataframe
    df = pyopenms.peptide_identifications_to_df(peps)
    df.loc[0, "ScoreType"] = 10.0
    peps = pyopenms.update_scores_from_df(peps, df, "ScoreType")
    assert peps[0].getHits()[0].getScore() == 10.0


if __name__ == "__main__":
    # Run all tests in this module
    import sys
    module = sys.modules[__name__]
    
    test_functions = [getattr(module, name) for name in dir(module) 
                     if name.startswith('test') and callable(getattr(module, name))]
    
    print(f"Running {len(test_functions)} identification tests...")
    for test_func in test_functions:
        try:
            test_func()
            print(f"✅ {test_func.__name__}")
        except Exception as e:
            print(f"❌ {test_func.__name__}: {e}")
