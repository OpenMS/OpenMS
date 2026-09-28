"""
Call forms that pyOpenMS 3.5.0 accepted and the nanobind port lost (#10260).

Each test calls the 3.5.0 form: trailing arguments that the C++ API defaults
are left out, and SimpleSearchEngineAlgorithm.search() is expected to return
its result instead of raising std::bad_cast.
"""

import pytest

import pyopenms as oms
from pyopenms.Constants import PROTON_MASS_U


def _write_search_input(tmp_path, peptide="DLGEEHFK"):
    """Write a FASTA containing peptide and an mzML with its theoretical MS2 spectrum."""
    fasta = tmp_path / "db.fasta"
    fasta.write_text(">P1\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFK" + peptide + "GLVLIAFSQYLQQCPFDEHVK\n")

    sequence = oms.AASequence.fromString(peptide)
    spectrum = oms.MSSpectrum()
    tsg = oms.TheoreticalSpectrumGenerator()
    tsg.getSpectrum(spectrum, sequence, 1, 1)
    spectrum.setMSLevel(2)
    spectrum.setRT(10.0)
    spectrum.setNativeID("scan=1")
    precursor = oms.Precursor()
    precursor.setCharge(2)
    precursor.setMZ((sequence.getMonoWeight() + 2 * PROTON_MASS_U) / 2)
    spectrum.setPrecursors([precursor])
    experiment = oms.MSExperiment()
    experiment.addSpectrum(spectrum)
    mzml = tmp_path / "spectra.mzML"
    oms.MzMLFile().store(str(mzml), experiment)
    return str(mzml), str(fasta)


def test_simple_search_engine_returns_exit_code_and_proteins(tmp_path):
    mzml, fasta = _write_search_input(tmp_path)
    peptide_ids = oms.PeptideIdentificationList()
    exit_code, protein_ids = oms.SimpleSearchEngineAlgorithm().search(mzml, fasta, peptide_ids)

    assert exit_code == oms.SimpleSearchEngineAlgorithm.ExitCodes.EXECUTION_OK
    assert len(protein_ids) == 1
    hits = [hit.getSequence().toString() for pid in peptide_ids for hit in pid.getHits()]
    assert "DLGEEHFK" in hits


def test_simple_search_engine_exit_codes_are_arithmetic():
    assert int(oms.SimpleSearchEngineAlgorithm.ExitCodes.EXECUTION_OK) == 0


def test_residue_weights_default_to_the_full_residue():
    residue = oms.ResidueDB().getResidue("A")
    full = oms.Residue.ResidueType.Full
    assert residue.getMonoWeight() == residue.getMonoWeight(full)
    assert residue.getAverageWeight() == residue.getAverageWeight(full)
    assert residue.getFormula() == residue.getFormula(full)


def test_residue_constructor_without_optional_properties():
    residue = oms.Residue("Alanine", "Ala", "A", oms.EmpiricalFormula("C3H7NO2"))
    assert residue.getOneLetterCode() == "A"


def test_fine_isotope_pattern_generator_takes_the_threshold_alone():
    formula = oms.EmpiricalFormula("C100H202")
    one_arg = oms.FineIsotopePatternGenerator(0.01).run(formula)
    two_args = oms.FineIsotopePatternGenerator(0.01, True).run(formula)
    full = oms.FineIsotopePatternGenerator(0.01, True, False).run(formula)
    assert one_arg.size() == two_args.size() == full.size() > 0


def test_meta_info_get_value_without_default():
    info = oms.MetaInfo()
    info.setValue("x", 5)
    assert info.getValue("x") == 5
    assert info.getValue("missing") is None


def test_modification_definition_with_mod_and_fixed():
    definition = oms.ModificationDefinition("Oxidation (M)", False)
    assert not definition.isFixedModification()
    assert oms.ModificationDefinition("Oxidation (M)").isFixedModification()


def test_id_filter_and_conflict_resolver_defaults():
    oms.IDFilter.removeDuplicatePeptideHits(oms.PeptideIdentificationList())
    oms.IDConflictResolverAlgorithm.resolve(oms.FeatureMap())
    oms.IDConflictResolverAlgorithm.resolve(oms.ConsensusMap())


def test_text_file_three_argument_form(tmp_path):
    path = tmp_path / "lines.txt"
    path.write_text("  a  \nb\n")
    text = oms.TextFile(str(path), True, -1)  # trim_lines=True
    out = tmp_path / "stored.txt"
    text.store(str(out))
    assert out.read_text().splitlines() == ["a", "b"]


def test_constructors_with_defaulted_arguments(tmp_path):
    oms.MzMLSpectrumDecoder()
    oms.IMSIsotopeDistribution()
    oms.BayesianProteinInferenceAlgorithm()
    oms.MSDataCachedConsumer(str(tmp_path / "data.cached"))
    oms.Tagger(2, 10.0, 65535, 1, 1, [], [])
    oms.DigestionEnzymeProtein("Test", "(?<=K)", set(), "", oms.EmpiricalFormula("H"),
                               oms.EmpiricalFormula("OH"), "", "", -1, -1)


def test_mass_explainer_compute_without_argument():
    explainer = oms.MassExplainer()
    explainer.compute()


def test_java_info_can_run_without_verbosity_flag():
    assert oms.JavaInfo.canRun("a-java-executable-that-does-not-exist") is False


def _identified_run(shift):
    """A FeatureMap and its PeptideIdentificationList: five peptides, RTs shifted by 'shift'."""
    fmap = oms.FeatureMap()
    peptides = oms.PeptideIdentificationList()
    for i, sequence in enumerate(["PEPTIDEK", "ELVISLIVESK", "DLGEEHFK", "LVNELTEFAK", "YLYEIAR"]):
        rt = 100.0 * (i + 1) + shift
        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(sequence))
        hit.setScore(0.01)
        pid = oms.PeptideIdentification()
        pid.setRT(rt)
        pid.setMZ(500.0 + i)
        pid.setScoreType("q-value")
        pid.setHigherScoreBetter(False)
        pid.setHits([hit])
        peptides.push_back(pid)
        feature = oms.Feature()
        feature.setRT(rt)
        feature.setMZ(500.0 + i)
        feature_peptides = oms.PeptideIdentificationList()
        feature_peptides.push_back(pid)
        feature.setPeptideIdentifications(feature_peptides)
        fmap.push_back(feature)
    return fmap, peptides


def _shifts(trafo):
    return sorted(round(p.first - p.second, 6) for p in trafo.getDataPoints())


def test_map_alignment_identification_aligns_several_maps():
    (fmap0, peps0), (fmap1, peps1) = _identified_run(0.0), _identified_run(10.0)
    aligner = oms.MapAlignmentAlgorithmIdentification()

    trafos = aligner.align([fmap0, fmap1], 0)
    assert len(trafos) == 2
    assert _shifts(trafos[1]) == [10.0] * 5

    # pyOpenMS 3.5 form: the transformations are written into the list passed in
    filled = [oms.TransformationDescription()]
    assert aligner.align([fmap0, fmap1], filled, 0) is None
    assert len(filled) == 2 and _shifts(filled[1]) == [10.0] * 5

    trafos = aligner.align([peps0, peps1], 0)
    assert _shifts(trafos[1]) == [10.0] * 5

    aligner.setReference(peps0)
    trafos = aligner.align([fmap1])
    assert len(trafos) == 1 and _shifts(trafos[0]) == [10.0] * 5


def test_id_filter_filters_a_protein_list_by_score():
    run = oms.ProteinIdentification()
    run.setHigherScoreBetter(True)
    hits = []
    for accession, score in (("P1", 0.9), ("P2", 0.5), ("P3", 0.1)):
        hit = oms.ProteinHit()
        hit.setAccession(accession)
        hit.setScore(score)
        hits.append(hit)
    run.setHits(hits)
    proteins = [run]
    filtered = oms.IDFilter.filterHitsByScore(proteins, 0.5)
    assert [h.getAccession() for h in filtered[0].getHits()] == ["P1", "P2"]
    assert [h.getAccession() for h in proteins[0].getHits()] == ["P1", "P2"]


def test_id_filter_update_protein_references_is_a_deprecated_alias():
    run = oms.ProteinIdentification()
    run.setIdentifier("run1")
    known = oms.ProteinHit()
    known.setAccession("P1")
    run.setHits([known])
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDEK"))
    evidences = []
    for accession in ("P1", "P2"):
        evidence = oms.PeptideEvidence()
        evidence.setProteinAccession(accession)
        evidences.append(evidence)
    hit.setPeptideEvidences(evidences)
    pid = oms.PeptideIdentification()
    pid.setIdentifier("run1")
    pid.setHits([hit])
    peptides = oms.PeptideIdentificationList()
    peptides.push_back(pid)
    with pytest.warns(DeprecationWarning, match="removeDanglingProteinReferences"):
        oms.IDFilter.updateProteinReferences(peptides, [run], False)
    remaining = [e.getProteinAccession() for e in peptides[0].getHits()[0].getPeptideEvidences()]
    assert remaining == ["P1"]


def test_opxl_helper_compute_delta_scores_takes_the_identifications():
    pid = oms.PeptideIdentification()
    pid.setHigherScoreBetter(True)
    hits = []
    for score in (5.0, 10.0):
        hit = oms.PeptideHit()
        hit.setScore(score)
        hits.append(hit)
    pid.setHits(hits)
    peptides = oms.PeptideIdentificationList()
    peptides.push_back(pid)
    assert oms.OPXLHelper.computeDeltaScores(peptides) is None
    scored = peptides[0].getHits()
    assert [h.getScore() for h in scored] == [10.0, 5.0]
    assert [h.getMetaValue("delta_score") for h in scored] == [0.5, 0.0]


def test_transition_files_accept_bytes_paths(tmp_path):
    target = tmp_path / "transitions.tsv"
    oms.TransitionTSVFile().convertTargetedExperimentToTSV(str(target).encode(), oms.TargetedExperiment())
    assert target.exists()
