"""
Call forms that pyOpenMS 3.5.0 accepted and the nanobind port lost (#10260).

Each test calls the 3.5.0 form: trailing arguments that the C++ API defaults
are left out, and SimpleSearchEngineAlgorithm.search() is expected to return
its result instead of raising std::bad_cast.
"""

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
