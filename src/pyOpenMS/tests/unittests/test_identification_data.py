"""Owning identification values, callback lifetimes and native streaming contracts."""

import copy
import gc
from pathlib import Path

import pytest
import pyopenms as oms


ID = oms.IdentificationData
File = oms.IdentificationDataFile


def make_run(name="search"):
    run = ID.Run(name)
    definition = ID.ScoreDefinition()
    definition.name = "posterior error probability"
    definition.higher_better = False
    definition.calibration = "test calibration"
    score = run.add_score(definition)
    source = ID.SourceFile()
    source.path = "/measurements/sample.mzML"
    source.primary_files = [source.path]
    source_id = run.add_source(source)
    observation = ID.Observation()
    observation.data_id = "scan=1"
    observation.rt = 12.5
    observation.mz = 400.2
    query = run.add_identification(source_id, observation)
    payload = ID.MatchData()
    payload.representation = "PEPTIDE"
    payload.charge = 2
    parent = ID.QualifiedAccession()
    parent.database = "db.fasta"
    parent.accession = "P1"
    evidence = ID.ParentEvidence()
    evidence.parent = parent
    evidence.start = 2
    evidence.end = 8
    payload.parent_evidence = [evidence]
    payload.setMetaValues({"empty": None})
    payload.setMetaValue("annotation", "")
    payload.setMetaValue("counts", [1, 2, 3])
    first = run.add_match(query, payload, [0.01])
    payload.representation = "EDITPEP"
    second = run.add_match(query, payload, [0.1])
    run.set_primary_score(score)
    run.set_selected_match(query, second)
    return run, query, first, second, score


def make_data_with_inference():
    run, query, first, second, score = make_run()
    data = ID()
    data.add_run(run)
    result = ID.InferenceResult()
    result.identifier = "pooled inference"
    protein = oms.ProteinHit()
    protein.setAccession("P1")
    protein.setScore(0.99)
    proteins = oms.ProteinIdentification()
    proteins.setScoreType("probability")
    proteins.setHigherScoreBetter(True)
    proteins.setHits([protein])
    result.proteins = proteins
    identity = ID.QualifiedAccession()
    identity.database = "db.fasta"
    identity.accession = "P1"
    result.qualified_accessions = {"P1": identity}
    input_record = ID.InferenceInput()
    input_record.run_identifier = run.get_identifier()
    input_record.run_uuid = run.get_uuid()
    input_record.score = run.get_score_definition(score)
    input_record.selection = "all candidates"
    result.inputs = [input_record]
    data.add_inference_result(result)
    return data, query, first, second


def test_score_contract_and_owned_records():
    run, query, first, _, score = make_run()
    retained = run.get_match(first)
    definition = run.get_score_definition(score)
    definition.name = "changed outside the owner"
    assert run.get_score_definition(score).name == "posterior error probability"
    retained.representation = "OTHER"
    assert run.get_match(first).representation == "PEPTIDE"
    with pytest.raises(Exception):
        run.set_score(first, score, None)
    with pytest.raises(Exception):
        run.set_score(first, score, float("nan"))
    assert run.get_score(first, score) == 0.01
    run.replace_match(first, retained, [0.01])
    assert run.get_match(first).representation == "OTHER"
    assert run.get_identification(query).get_matches()[0].get_id() == first
    assert not hasattr(run, "get_revision")
    assert not hasattr(ID(), "is_current")


def test_get_run_copy_can_be_committed_without_dangling_references():
    run, _, first, _, _ = make_run()
    data = ID()
    data.add_run(run)
    retained = data.get_run("search")
    for index in range(30):
        data.add_run(ID.Run(f"other-{index}"))
    retained.erase_matches(lambda match: match.get_id() == first)
    assert data.get_run("search").get_number_of_matches() == 2
    data.replace_run(retained)
    assert data.get_run("search").get_number_of_matches() == 1
    assert run.get_number_of_matches() == 2
    del data
    gc.collect()
    assert retained.get_number_of_matches() == 1


def test_nested_field_values_remain_owned_after_replacement():
    result = ID.InferenceResult()
    first_input = ID.InferenceInput()
    first_input.run_identifier = "first"
    result.inputs = [first_input]
    retained = result.inputs[0]
    result.inputs = []
    assert retained.run_identifier == "first"
    run, _, first, _, _ = make_run()
    evidence = run.get_match(first).parent_evidence[0]
    parents = run.get_match(first).parent_evidence
    parents.clear()
    assert evidence.parent.accession == "P1"
    assert len(run.get_match(first).parent_evidence) == 1


def test_copy_preserves_ids_and_edits_are_independent():
    run, query, first, _, score = make_run()
    clone = copy.deepcopy(run)
    assert clone.get_uuid() == run.get_uuid()
    assert clone.get_match(first).get_id() == first
    clone.set_score(first, clone.get_score_id(score.value), 0.2)
    assert run.get_score(first, score) == 0.01
    assert clone.get_score(first, clone.get_primary_score()) == 0.2
    assert clone.get_identification(query).data_id == "scan=1"


def test_callbacks_retain_safe_values_and_rollback_on_exceptions():
    run, _, first, _, score = make_run()
    retained = []

    def invalid_callback(match):
        retained.append(match)
        run.set_score(first, score, 0.5)
        return True

    with pytest.raises(Exception):
        run.filter_matches(invalid_callback)
    assert run.get_score(first, score) == 0.01
    assert run.get_number_of_matches() == 2

    def throwing_callback(match):
        retained.append(match)
        if len(retained) > 2:
            raise ValueError("callback failed")
        return False

    with pytest.raises(ValueError, match="callback failed"):
        run.filter_matches(throwing_callback)
    assert run.get_number_of_matches() == 2
    run.erase_matches(lambda match: True)
    del run
    gc.collect()
    assert retained[0].representation == "PEPTIDE"
    assert retained[-1].get_scores() == [0.1]


def test_transform_callbacks_use_owned_payloads_and_commit():
    run, _, first, _, _ = make_run()
    retained = []

    def transform(payload):
        retained.append(payload)
        payload.name = "updated"

    run.transform_matches(transform)
    assert run.get_match(first).name == "updated"
    retained[0].name = "changed later"
    assert run.get_match(first).name == "updated"
    run.erase_matches(lambda match: True)
    assert retained[1].representation == "EDITPEP"


def test_filter_preserves_inference_provenance_or_explicitly_discards():
    data, query, first, second = make_data_with_inference()
    data.filter_matches(lambda match: match.get_id() == first, ID.InferencePolicy.PRESERVE)
    run = data.get_run("search")
    assert run.find_match(second) is None
    assert run.get_identification(query).get_selected_match() is None
    result = data.get_inference_results()[0]
    assert result.inputs[0].run_uuid == run.get_uuid()
    assert result.inputs[0].selection == "all candidates"
    assert result.proteins.getHits()[0].getScore() == 0.99
    data.filter_matches(lambda match: True, ID.InferencePolicy.DISCARD)
    assert data.get_inference_results() == []


def test_native_round_trip_and_nonreused_ids(tmp_path):
    data, query, first, second = make_data_with_inference()
    data.filter_matches(lambda match: match.get_id() == first, ID.InferencePolicy.PRESERVE)
    path = str(tmp_path / "native")
    options = File.Options()
    options.batch_rows = 1
    options.row_group_rows = 1
    File.store(path, data, options)
    loaded = File.load(path)
    run = loaded.get_run("search")
    assert run.get_uuid() == data.get_run("search").get_uuid()
    match = run.get_match(first)
    assert match.representation == "PEPTIDE"
    assert match.metaValueExists("empty")
    assert match.getMetaValue("empty") is None
    assert match.getMetaValue("annotation") == ""
    assert match.getMetaValue("counts") == [1, 2, 3]
    assert match.parent_evidence[0].parent.database == "db.fasta"
    result = loaded.get_inference_results()[0]
    assert result.inputs[0].run_uuid == run.get_uuid()
    assert result.inputs[0].selection == "all candidates"
    added = run.add_match(query, match, [0.05])
    assert added.value > second.value
    descriptor = File.inspect(path)[0]
    assert descriptor.match_count == 1
    assert descriptor.query_count == 1
    assert File.load_run(path, descriptor.uuid).get_number_of_matches() == 1


def test_native_projected_scan_keeps_callback_values_alive(tmp_path):
    run, _, _, _, _ = make_run()
    data = ID()
    data.add_run(run)
    path = str(tmp_path / "native")
    options = File.Options()
    options.batch_rows = 1
    options.row_group_rows = 1
    File.store(path, data, options)
    scan = File.ScanOptions()
    scan.buffering = options
    scan.runs = [run.get_uuid()]
    projection = scan.projection
    projection.molecule = False
    projection.evidence = False
    projection.annotations = False
    projection.metadata = False
    projection.all_scores = False
    projection.score_ids = [0]
    scan.projection = projection
    batches = []
    queries = []
    stats = File.scan(path, scan, lambda uuid, batch: queries.extend(batch), lambda uuid, batch: batches.append(batch))
    assert stats.matches == 2
    assert stats.queries == 1
    assert stats.descriptor_bytes > 0
    assert len(batches) == 2
    assert [batch[0].scores for batch in batches] == [[0.01], [0.1]]
    assert batches[0][0].data.representation == ""
    assert batches[0][0].data.parent_evidence == []
    assert queries[0].data.data_id == "scan=1"
    gc.collect()
    assert batches[0][0].match_id == 1


def test_native_filter_writes_reduced_dataset_with_explicit_inference_policy(tmp_path):
    data, query, first, second = make_data_with_inference()
    source = str(tmp_path / "source")
    destination = str(tmp_path / "reduced")
    File.store(source, data)
    File.filter(source, destination, lambda uuid, record: record.scores[0] < 0.05, ID.InferencePolicy.PRESERVE)
    result = File.load(destination)
    assert result.get_run("search").get_number_of_matches() == 1
    assert result.get_run("search").get_match(first).get_id() == first
    assert result.get_run("search").get_identification(query).get_selected_match() is None
    assert result.get_inference_results()[0].inputs[0].run_uuid == data.get_run("search").get_uuid()
    assert result.get_inference_results()[0].inputs[0].selection == "all candidates"
    stripped = str(tmp_path / "stripped")
    File.filter(destination, stripped, lambda uuid, record: True, ID.InferencePolicy.DISCARD)
    assert File.load(stripped).get_inference_results() == []
    with pytest.raises(Exception):
        File.store(source, data)


def test_failed_native_load_preserves_destination(tmp_path):
    run, _, first, _, _ = make_run()
    data = ID()
    data.add_run(run)
    path = str(tmp_path / "native")
    File.store(path, data)
    matches = list(Path(path).glob("matches.parquet"))
    assert len(matches) == 1
    matches[0].unlink()
    with pytest.raises(Exception):
        File.load(path, data)
    assert data.get_run("search").get_match(first).representation == "PEPTIDE"
    assert data.get_run("search").get_uuid() == run.get_uuid()


def test_optional_adduct_and_empty_queries_round_trip(tmp_path):
    run = ID.Run("compounds", ID.MoleculeKind.COMPOUND)
    score = ID.ScoreDefinition()
    score.name = "similarity"
    run.set_primary_score(run.add_score(score))
    source = run.add_source(ID.SourceFile())
    query = run.add_identification(source, ID.Observation())
    empty = run.add_identification(source, ID.Observation())
    payload = ID.MatchData()
    payload.encoding = ID.Encoding.SMILES
    payload.representation = "CCO"
    payload.charge = 1
    payload.adduct = oms.AMSE_AdductInfo.parseAdductString("2M+H;1+")
    match_id = run.add_match(query, payload, [0.9])
    data = ID()
    data.add_run(run)
    path = str(tmp_path / "native")
    File.store(path, data)
    loaded = File.load(path).get_run("compounds")
    assert loaded.get_number_of_identifications() == 2
    assert loaded.get_identification(empty).get_matches() == []
    adduct = loaded.get_match(match_id).adduct
    assert adduct.getCharge() == 1
    assert adduct.getMolMultiplier() == 2
    assert adduct.getName() == "2M+H;1+"


def test_legacy_adapter_preserves_scores_and_molecular_values():
    protein = oms.ProteinHit()
    protein.setAccession("P1")
    proteins = oms.ProteinIdentification()
    proteins.setIdentifier("legacy-search")
    proteins.setHits([protein])
    evidence = oms.PeptideEvidence()
    evidence.setProteinAccession("P1")
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setScore(0.03)
    hit.setPeptideEvidences([evidence])
    query = oms.PeptideIdentification()
    query.setIdentifier("legacy-search")
    query.setScoreType("PEP")
    query.setHigherScoreBetter(False)
    query.setHits([hit])
    peptides = oms.PeptideIdentificationList()
    peptides.append(query)
    imported = oms.IdentificationDataAdapter.import_legacy([proteins], peptides)
    assert len(imported.queries) == 1
    values = imported.data
    assert values.get_runs()[0].get_number_of_matches() == 1
    exported = oms.IdentificationDataAdapter.to_legacy(values)
    assert exported.losses == []
    assert exported.peptides[0].getHits()[0].getSequence().toString() == "PEPTIDE"
    assert exported.peptides[0].getHits()[0].getScore() == 0.03
    assert len(exported.queries) == 1


def test_inference_pools_runs_without_editing_candidates():
    data = ID()
    inputs = []
    for name in ("analysis-A", "analysis-B"):
        run, _, _, _, score = make_run(name)
        data.add_run(run)
        item = oms.IdentificationDataInference.Input()
        item.run_uuid = run.get_uuid()
        item.score = score
        inputs.append(item)
    result = oms.IdentificationDataInference.infer(data, inputs, "pooled")
    assert len(result.inputs) == 2
    assert [entry.run_uuid for entry in result.inputs] == [entry.run_uuid for entry in inputs]
    assert len(result.proteins.getHits()) == 1
    assert [run.get_number_of_matches() for run in data.get_runs()] == [2, 2]
    assert result.parent_score.scope == ID.ScoreScope.PROTEIN
    assert data.get_inference_results() == []
    data.add_inference_result(result)
    oms.IdentificationDataInference.retain_proteins(result, [])
    assert result.proteins.getHits() == []
    assert len(data.get_inference_results()[0].proteins.getHits()) == 1
    assert data.get_run("analysis-A").get_number_of_matches() == 2


def test_file_handler_loads_and_stores_native_owning_values(tmp_path):
    run, _, first, _, _ = make_run()
    data = ID()
    data.add_run(run)
    path = str(tmp_path / "search.idparquet")
    handler = oms.FileHandler()
    handler.store_identification_data(path, data)
    loaded = handler.load_identification_data(path)
    assert loaded.get_run("search").get_uuid() == run.get_uuid()
    assert loaded.get_run("search").get_match(first).representation == "PEPTIDE"


def test_dataset_ordered_score_contract():
    first, _, _, _, _ = make_run("A")
    second, _, _, _, _ = make_run("B")
    extra = ID.ScoreDefinition()
    extra.name = "supplementary"
    first.add_score(extra)
    second.add_score(extra)
    data = ID()
    data.add_run(first)
    data.add_run(second)
    definition = data.get_primary_score_definition()
    data.set_primary_score(definition)
    data.validate()
    assert [score.name for score in data.get_score_definitions()] == [definition.name, "supplementary"]
    with pytest.raises(Exception):
        data.set_primary_score(extra)  # Null supplementary values cannot become primary.
    assert data.get_primary_score_definition().name == definition.name
    reordered = ID.Run("reordered")
    reordered.add_score(extra)
    reordered.set_primary_score(reordered.add_score(definition))
    with pytest.raises(Exception):
        data.add_run(reordered)
    incompatible, _, _, _, _ = make_run("bad")
    other = ID.ScoreDefinition()
    other.name = "different score"
    score_id = incompatible.add_score(other)
    for block in incompatible.get_source_blocks():
        for query in block.identifications:
            for match in query.get_matches():
                incompatible.set_score(match.get_id(), score_id, 1.0)
    incompatible.set_primary_score(score_id)
    with pytest.raises(Exception):
        data.add_run(incompatible)
    assert len(data.get_runs()) == 2
    with pytest.raises(Exception):
        data.set_primary_score(other)
    assert data.get_primary_score_definition().name == definition.name
