"""In/out list arguments update the caller's list, as in pyOpenMS 3.5, and are also returned.

apply_addons() installs pyopenms/addons/inout_lists.py when pyopenms is imported.
"""
import re

import pytest

import pyopenms as oms
from pyopenms.addons import inout_lists


# --------------------------------------------------------------------- helpers
def _protein_run(accessions):
    run = oms.ProteinIdentification()
    run.setIdentifier("run1")
    hits = []
    for i, acc in enumerate(accessions):
        h = oms.ProteinHit()
        h.setAccession(acc)
        h.setScore(float(len(accessions) - i))
        h.setMetaValue("target_decoy", "decoy" if acc.startswith("DECOY_") else "target")
        hits.append(h)
    run.setHits(hits)
    run.setHigherScoreBetter(True)
    return run


def _peptides(accessions, identifier="run1"):
    peps = oms.PeptideIdentificationList()
    for acc in accessions:
        ev = oms.PeptideEvidence()
        ev.setProteinAccession(acc)
        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString("PEPTIDER"))
        hit.setScore(0.01)
        hit.setPeptideEvidences([ev])
        hit.setMetaValue("target_decoy", "decoy" if acc.startswith("DECOY_") else "target")
        pid = oms.PeptideIdentification()
        pid.setIdentifier(identifier)
        pid.setScoreType("q-value")
        pid.setHigherScoreBetter(False)
        pid.setHits([hit])
        peps.push_back(pid)
    return peps


def _accessions(prots):
    return [h.getAccession() for h in prots[0].getHits()]


# ---------------------------------------------------------- table consistency
_SIG = re.compile(r"def \w+\((?P<params>.*)\) -> (?P<ret>.+)$")


def _split_params(text):
    out, depth, cur = [], 0, ""
    for ch in text:
        depth += ch in "[("
        depth -= ch in "])"
        if ch == "," and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return [p for p in out if p not in ("self", "/", "*")]


def _is_sequence(annotation):
    # std::vector<T> renders as Sequence[T]; OpenMS' own std::vector<std::string> caster as list[str]
    return annotation.startswith(("collections.abc.Sequence[", "list["))


def _overloads(cls, meth):
    raw = cls.__dict__[meth]
    func = raw.__func__ if isinstance(raw, staticmethod) else raw
    for sig, _doc, _defaults in func.__nb_signature__:
        m = _SIG.match(sig.replace("\n", " "))
        params = _split_params(m.group("params"))
        yield [p.split(":")[0].strip() for p in params], [p.split(":", 1)[1].split("=")[0].strip() for p in params], m.group("ret").strip()


@pytest.mark.parametrize("cls_name,meth", sorted(inout_lists.INOUT_METHODS))
def test_table_matches_binding_signatures(cls_name, meth):
    """Every rule must name a list parameter of an overload that returns that shape."""
    rules = inout_lists.INOUT_METHODS[(cls_name, meth)]
    cls = getattr(oms, cls_name)
    overloads = list(_overloads(cls, meth))
    if inout_lists.LIST in rules:
        pos, names = rules[inout_lists.LIST]
        assert any(r.startswith("list[") and len(n) > pos and n[pos] in names and _is_sequence(t[pos])
                   for n, t, r in overloads), overloads
    for length, mapping in rules.get(inout_lists.TUPLE, {}).items():
        for _idx, pos, names in mapping:
            assert any(r.startswith("tuple") and len(n) > pos and n[pos] in names and _is_sequence(t[pos])
                       for n, t, r in overloads), overloads


def test_install_is_idempotent_and_keeps_metadata():
    assert inout_lists.install(oms) == []
    f = oms.IDFilter.__dict__["removeUnreferencedProteins"]
    assert isinstance(f, staticmethod)
    assert oms.IDFilter.removeUnreferencedProteins.__name__ == "removeUnreferencedProteins"
    doc = oms.IDFilter.removeUnreferencedProteins.__doc__
    assert "removeUnreferencedProteins(cmap" in doc and "updated in place" in doc
    assert oms.FalseDiscoveryRate.apply.__qualname__ == "FalseDiscoveryRate.apply"


# ------------------------------------------------------------------ behaviour
@pytest.mark.parametrize("call", ["class", "instance", "keywords", "ids_keyword"])
def test_removeUnreferencedProteins_both_call_styles(call):
    prots = [_protein_run(["P1", "P2", "P3"])]
    peps = _peptides(["P2"])
    if call == "class":
        ret = oms.IDFilter.removeUnreferencedProteins(prots, peps)
    elif call == "instance":
        ret = oms.IDFilter().removeUnreferencedProteins(prots, peps)
    elif call == "keywords":
        ret = oms.IDFilter.removeUnreferencedProteins(proteins=prots, peptides=peps)
    else:
        ret = oms.IDFilter.removeUnreferencedProteins(prots, ids=peps)
    assert _accessions(ret) == ["P2"]      # 3.6 style: returned copy
    assert _accessions(prots) == ["P2"]    # 3.5 style: caller's list updated
    assert ret is not prots                # return value is still a separate list


def test_non_list_sequence_is_accepted_and_left_alone():
    prots = (_protein_run(["P1", "P2"]),)
    ret = oms.IDFilter.removeUnreferencedProteins(prots, _peptides(["P1"]))
    assert _accessions(ret) == ["P1"]
    assert _accessions(prots) == ["P1", "P2"]


def test_none_returning_overload_untouched():
    """FalseDiscoveryRate.apply(PeptideIdentificationList, bool) updates in place and returns None."""
    peps = _peptides(["P1", "DECOY_P2"])
    assert oms.FalseDiscoveryRate().apply(peps, False) is None
    assert peps[0].getScoreType() == "q-value"


def test_fdr_apply_protein_list_and_tuple_forms():
    prots = [_protein_run(["P1", "DECOY_P2", "P3"])]
    ret = oms.FalseDiscoveryRate().apply(prots)
    assert prots[0].getScoreType() == ret[0].getScoreType() == "q-value"
    assert _accessions(prots) == _accessions(ret)

    fwd, rev = [_protein_run(["P1", "P2"])], [_protein_run(["DECOY_P1"])]
    ret = oms.FalseDiscoveryRate().apply(fwd, rev)
    assert isinstance(ret, tuple) and len(ret) == 2
    assert fwd[0].getScoreType() == ret[0][0].getScoreType() == "q-value"
    assert [h.getScore() for h in fwd[0].getHits()] == [h.getScore() for h in ret[0][0].getHits()]


def test_basic_protein_inference_run():
    prots = [_protein_run(["P1", "P2"])]
    peps = _peptides(["P1", "P1", "P2"])
    ret = oms.BasicProteinInferenceAlgorithm().run(peps, prots)
    assert [h.getScore() for h in prots[0].getHits()] == [h.getScore() for h in ret[0].getHits()]
    assert prots[0].getScoreType() == ret[0].getScoreType() != ""


def test_transformation_model_weight_data():
    p = oms.Param()
    p.setValue("x_weight", "ln(x)")
    p.setValue("y_weight", "y")
    data = [oms.TM_DataPoint(10.0, 1.0), oms.TM_DataPoint(100.0, 2.0), oms.TM_DataPoint(1000.0, 3.0)]
    model = oms.TransformationModelLinear(data, p)
    ret = model.weightData(data)
    assert [round(d.first, 4) for d in data] == [round(d.first, 4) for d in ret] == [2.3026, 4.6052, 6.9078]
    model.unWeightData(data)
    assert [round(d.first, 4) for d in data] == [10.0, 100.0, 1000.0]


@pytest.mark.parametrize("keyword", [False, True])
def test_percolator_comet_features(keyword):
    fs = ["existing"]
    helper = oms.PercolatorFeatureSetHelper
    if keyword:
        ret = helper().addCOMETFeatures(oms.PeptideIdentificationList(), feature_set=fs)
    else:
        ret = helper.addCOMETFeatures(oms.PeptideIdentificationList(), fs)
    assert fs == ret and fs[0] == "existing" and len(fs) == 10


def _opxl_peptide_ids():
    pid = oms.PeptideIdentification()
    pid.setHigherScoreBetter(True)
    hits = []
    for score in (5.0, 10.0):
        hit = oms.PeptideHit()
        hit.setScore(score)
        hit.setMetaValue("target_decoy", "target")
        hits.append(hit)
    pid.setHits(hits)
    return pid


def test_opxl_peptide_identification_list_is_updated():
    pil = oms.PeptideIdentificationList()
    pil.push_back(_opxl_peptide_ids())
    assert oms.OPXLHelper.addXLTargetDecoyMV(pil) is None
    assert pil[0].getHits()[0].metaValueExists("xl_target_decoy_alpha")


def test_opxl_python_list_is_updated():
    peptide_ids = [_opxl_peptide_ids()]
    ret = oms.OPXLHelper.addXLTargetDecoyMV(peptide_ids)
    assert ret[0].getHits()[0].metaValueExists("xl_target_decoy_alpha")
    assert peptide_ids[0].getHits()[0].metaValueExists("xl_target_decoy_alpha")


def test_error_leaves_list_unchanged():
    prots = [_protein_run(["P1"])]
    with pytest.raises(TypeError):
        oms.IDFilter.removeUnreferencedProteins(prots, "not a peptide list")
    assert _accessions(prots) == ["P1"]


def test_pep_model_fit_probabilities_written_back():
    scores = [0.1 * i for i in range(-50, 50)] + [3.0 + 0.05 * i for i in range(40)]
    probs = []
    ok, sorted_scores, ret_probs = oms.PosteriorErrorProbabilityModel().fit(scores, probs, "none")
    assert probs == ret_probs and len(probs) == len(scores)
    assert scores == sorted_scores
