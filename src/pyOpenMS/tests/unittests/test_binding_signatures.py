"""Bound functions only use types that have a Python binding.

A parameter or return type without a binding shows up as its raw C++ name in
the nanobind signature ("OpenMS::ConsensusFeature::Ratio", "std::fpos<...>") and
makes the function uncallable: nanobind can neither convert such an argument
nor such a return value.
"""
import gc
import inspect
import os
import re

import pytest

import pyopenms as oms

_CPP_TYPE = re.compile(r"\b(OpenMS::|std::|__gnu_cxx::|boost::|Eigen::)")
_DATA = os.path.dirname(os.path.abspath(__file__))


def _classes(name, cls, seen):
    """cls and the classes nested in it ("ConsensusFeature.Ratio"), once each."""
    if id(cls) in seen:
        return
    seen.add(id(cls))
    yield name, cls
    for attr, value in vars(cls).items():
        # only classes defined here, not aliases of other classes
        if inspect.isclass(value) and value.__qualname__ == f"{cls.__qualname__}.{attr}":
            yield from _classes(f"{name}.{attr}", value, seen)


def _signatures():
    seen = set()
    for name in dir(oms):
        obj = getattr(oms, name)
        # the DB singletons (ElementDB, ...) are exported as callable wrappers
        obj = getattr(obj, "_singleton_class", obj)
        if not inspect.isclass(obj):
            for sig, _doc, _defaults in getattr(obj, "__nb_signature__", None) or ():
                yield name, sig
            continue
        for where, cls in _classes(name, obj, seen):
            for attr, value in vars(cls).items():
                funcs = [value.fget, value.fset] if isinstance(value, property) else [value]
                for func in funcs:
                    func = getattr(func, "__func__", func)
                    for sig, _doc, _defaults in getattr(func, "__nb_signature__", None) or ():
                        yield f"{where}.{attr}", sig


def test_no_raw_cpp_types_in_signatures():
    offenders = sorted({f"{where}: {sig}" for where, sig in _signatures() if _CPP_TYPE.search(sig)})
    assert not offenders, "Signatures with unbound C++ types:\n" + "\n".join(offenders)


def test_stream_offsets_are_ints(tmp_path):
    offset = oms.IndexedMzMLDecoder().findIndexListOffset(os.path.join(_DATA, "test.indexed.mzML"), 1024)
    assert isinstance(offset, int) and offset > 0

    exp = oms.MSExperiment()
    for rt in (1.0, 2.0):
        spec = oms.MSSpectrum()
        spec.setRT(rt)
        spec.set_peaks(([100.0, 200.0], [1.0, 2.0]))
        exp.addSpectrum(spec)
    cached = str(tmp_path / "test.cached.mzML")
    handler = oms.CachedMzMLHandler()
    handler.writeMemdump(exp, cached)
    handler.createMemdumpIndex(cached)
    index = handler.getSpectraIndex()
    assert len(index) == 2 and all(isinstance(i, int) for i in index)
    assert handler.getChromatogramIndex() == []


def test_consensus_feature_ratios():
    ratio = oms.ConsensusFeature.Ratio()
    assert ratio.ratio_value_ == 0.0
    ratio.ratio_value_ = 2.5
    ratio.numerator_ref_ = "heavy"
    ratio.denominator_ref_ = "light"
    ratio.description_ = ["SILAC"]
    feature = oms.ConsensusFeature()
    feature.addRatio(ratio)
    feature.setRatios(feature.getRatios() + [oms.ConsensusFeature.Ratio()])
    ratios = feature.getRatios()
    assert len(ratios) == 2
    assert (ratios[0].ratio_value_, ratios[0].numerator_ref_, ratios[0].description_) == (2.5, "heavy", ["SILAC"])


def test_mz_trafo_model_coefficients():
    model = oms.MZTrafoModel()
    with pytest.raises(Exception):
        model.getCoefficients()  # not trained yet
    model.setCoefficients(1.0, 2.0, 3.0)
    assert model.getCoefficients() == (1.0, 2.0, 3.0)


def test_noise_estimator():
    mz = [100.0 + i for i in range(100)]
    intensity = [10.0 + (i % 7) for i in range(100)]
    noise = oms.SignalToNoiseEstimatorMedianRapid(20.0).estimateNoise(mz, intensity)
    assert isinstance(noise, oms.SignalToNoiseEstimatorMedianRapid.NoiseEstimator)
    assert noise.get_noise_value(150.0) >= 1.0
    with pytest.raises(IndexError):
        noise.get_noise_value(10000.0)
    with pytest.raises(IndexError):
        noise.get_noise_even(float("nan"))
    # each accessor checks only the window set it reads
    estimator = oms.SignalToNoiseEstimatorMedianRapid.NoiseEstimator(3, 0, 10)
    estimator.result_windows_even = [4.0, 5.0, 6.0]
    estimator.result_windows_odd = []
    assert estimator.get_noise_even(5.0) == 4.0
    with pytest.raises(IndexError):
        estimator.get_noise_odd(5.0)
    with pytest.raises(IndexError):
        estimator.get_noise_value(5.0)


def test_spline_navigator_keeps_spline_alive():
    mz = [100.0 + 0.01 * i for i in range(50)]
    intensity = [1000.0 - (i - 25) ** 2 for i in range(50)]
    navigator = oms.SplineInterpolatedPeaks(mz, intensity).getNavigator()
    gc.collect()
    assert navigator.eval(100.25) == pytest.approx(1000.0, rel=1e-2)
    assert navigator.getNextPos(100.25) > 100.25


def test_isotope_labeling_mdvs_enums():
    feature = oms.Feature()
    subordinates = []
    for intensity in (1.0, 2.0, 4.0):
        sub = oms.Feature()
        sub.setMetaValue("peak_apex_int", intensity)
        subordinates.append(sub)
    feature.setSubordinates(subordinates)
    mdvs = oms.IsotopeLabelingMDVs()
    norm_max = mdvs.calculateMDV(feature, oms.IsotopeLabelingMDVs.MassIntensityType.NORM_MAX, "peak_apex_int")
    norm_sum = mdvs.calculateMDV(feature, oms.IsotopeLabelingMDVs.MassIntensityType.NORM_SUM, "peak_apex_int")
    assert [s.getMetaValue("peak_apex_int") for s in norm_max.getSubordinates()] == pytest.approx([0.25, 0.5, 1.0])
    assert [s.getMetaValue("peak_apex_int") for s in norm_sum.getSubordinates()] == pytest.approx([1 / 7, 2 / 7, 4 / 7])
    assert oms.IsotopeLabelingMDVs.DerivatizationAgent.TBDMS != oms.IsotopeLabelingMDVs.DerivatizationAgent.NOT_SELECTED


def test_ims_weights():
    weights = oms.IMSWeights([1.0, 2.5], 0.5)
    assert len(weights) == 2
    assert (weights.getWeight(0), weights.getWeight(1)) == (2, 5)
    assert weights.getAlphabetMass(1) == 2.5
    assert weights.getParentMass([1, 2]) == pytest.approx(6.0)
    with pytest.raises(IndexError):
        weights.getWeight(2)
    weights.swap(0, 1)
    assert weights.getAlphabetMass(0) == 2.5
    oms.RealMassDecomposer(weights)


def test_ims_weights_reject_what_the_cpp_divides_by():
    # each of these crashed the interpreter (division by zero, empty alphabet)
    with pytest.raises(ValueError):
        oms.IMSWeights([1.0, 2.0], 0.0)
    with pytest.raises(ValueError):
        oms.IMSWeights([1.0, 2.0], 0.5).setPrecision(-1.0)
    assert oms.IMSWeights([0.1, 0.2], 1.0).divideByGCD() is False  # all weights round to 0
    with pytest.raises(ValueError):
        oms.RealMassDecomposer(oms.IMSWeights())
    with pytest.raises(ValueError):
        oms.RealMassDecomposer(oms.IMSWeights([0.1, 57.02], 1.0))  # 0.1 rounds to weight 0


def test_picked_mass_trace():
    trace = oms.FeatureFinderAlgorithmPickedHelperStructs_MassTrace()
    trace.max_rt = 12.5
    assert oms.FeatureFinderAlgorithmPickedHelperStructs_MassTrace(trace).max_rt == 12.5
    assert not trace.isValid()


def test_delta_mass_label_set():
    delta = oms.MultiplexDeltaMasses_DeltaMass(8.0142, ["Lys8", "Arg10", "Lys8"])
    assert delta.label_set == ["Arg10", "Lys8", "Lys8"]
    delta.label_set = ["Arg6"]
    assert delta.label_set == ["Arg6"]
    # a str is one label, not a list of characters
    assert oms.MultiplexDeltaMasses_DeltaMass(8.0142, "Lys8").label_set == ["Lys8"]
    delta.label_set = "Arg10"
    assert delta.label_set == ["Arg10"]
    with pytest.raises(TypeError):
        delta.label_set = 5


def test_singleton_wrappers_stay_introspectable():
    # help() and IPython call inspect.signature() on ElementDB, ModificationsDB, ...
    assert str(inspect.signature(oms.ElementDB)) == "()"
    assert oms.ElementDB() is not None


def test_modification_map_size():
    mods = oms.ModifiedPeptideGenerator.getModifications(["Oxidation (M)", "Carbamidomethyl (C)"])
    assert len(mods) == 2


def test_peptide_indexing_from_fasta_file(tmp_path):
    fasta = tmp_path / "db.fasta"
    fasta.write_text(">P1\nMKPEPTIDERK\n>P2\nAAAAAAK\n")
    prot = oms.ProteinIdentification()
    prot.setIdentifier("run")
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDER"))
    pep = oms.PeptideIdentification()
    pep.setIdentifier("run")
    pep.setHits([hit])
    peps = oms.PeptideIdentificationList()
    peps.push_back(pep)
    indexer = oms.PeptideIndexing()
    params = indexer.getParameters()
    params.setValue("decoy_string", "DECOY_")
    params.setValue("missing_decoy_action", "warn")
    indexer.setParameters(params)
    _exit_code, prot_ids = indexer.run(str(fasta), [prot], peps)
    assert [h.getAccession() for h in prot_ids[0].getHits()] == ["P1"]
    assert [e.getProteinAccession() for e in peps[0].getHits()[0].getPeptideEvidences()] == ["P1"]
