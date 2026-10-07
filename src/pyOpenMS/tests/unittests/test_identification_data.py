"""Owning identification values, callback lifetimes and native streaming contracts."""

import copy
import gc
from pathlib import Path
import subprocess
import sys

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
    score = run.addScore(definition)
    source = ID.SourceFile()
    source.path = "/measurements/sample.mzML"
    source_id = run.addSource(source)
    observation = ID.Observation()
    observation.data_id = "scan=1"
    observation.rt = 12.5
    observation.mz = 400.2
    query = run.addIdentification(source_id, observation)
    payload = ID.MatchData()
    payload.representation = "PEPTIDE"
    payload.charge = 2
    database = ID.Database()
    database.path = "db.fasta"
    evidence = ID.SequenceEvidence()
    evidence.database = run.addDatabase(database)
    evidence.accession = "P1"
    evidence.start = 2
    evidence.end = 8
    payload.sequence_evidence = [evidence]
    payload.setMetaValues({"empty": None})
    payload.setMetaValue("annotation", "")
    payload.setMetaValue("counts", [1, 2, 3])
    first = run.addMatch(query, payload, [0.01])
    payload.representation = "EDITPEP"
    second = run.addMatch(query, payload, [0.1])
    run.setPrimaryScore(score)
    run.setSelectedMatch(query, second)
    return run, query, first, second, score


def make_data_with_inference():
    run, query, first, second, score = make_run()
    data = ID()
    data.addRun(run)
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
    input_record.run_identifier = run.getIdentifier()
    input_record.run_uuid = run.getUuid()
    input_record.score = run.getScoreDefinition(score)
    input_record.selection = "all candidates"
    result.inputs = [input_record]
    data.addInferenceResult(result)
    return data, query, first, second


def test_score_contract_and_owned_records():
    run, query, first, _, score = make_run()
    retained = run.getMatch(first)
    definition = run.getScoreDefinition(score)
    definition.name = "changed outside the owner"
    assert run.getScoreDefinition(score).name == "posterior error probability"
    retained.representation = "OTHER"
    assert run.getMatch(first).representation == "PEPTIDE"
    with pytest.raises(Exception):
        run.setScore(first, score, None)
    with pytest.raises(Exception):
        run.setScore(first, score, float("nan"))
    assert run.getScore(first, score) == 0.01
    run.replaceMatch(first, retained, [0.01])
    assert run.getMatch(first).representation == "OTHER"
    assert run.getIdentification(query).getMatches()[0].getId() == first
    assert not hasattr(run, "get_revision")
    assert not hasattr(ID(), "is_current")


def test_run_views_edit_the_live_run_and_copies_stay_snapshots():
    run, query, first, second, _ = make_run()
    data = ID()
    view = data.addRun(run)
    assert isinstance(view, ID.RunView)
    snapshot = data.getRun("search")
    for index in range(30):
        data.addRun(ID.Run(f"other-{index}"))  # growing the dataset does not invalidate the view
    view.eraseMatches(lambda match: match.getId() == first)
    assert data.getRun("search").getNumberOfMatches() == 1
    assert snapshot.getNumberOfMatches() == 2  # a copy is a snapshot; nothing writes it back
    assert not hasattr(data, "replaceRun") and not hasattr(data, "replace_run")
    # New IDs come from the live run only, so they never collide with IDs handed out before.
    added = data.run_view("search").addMatch(query, snapshot.getMatch(first), [0.3])
    assert added.value == second.value + 1
    assert view.getMatch(added).representation == "PEPTIDE"
    assert data.run_view_by_uuid(run.getUuid()).getNumberOfMatches() == 2
    assert view.copy().getNumberOfMatches() == 2
    # The view keeps the dataset alive and reports a run that is gone.
    del data
    gc.collect()
    assert view.getNumberOfMatches() == 2
    owner = ID()
    gone = owner.addRun("temporary")
    owner.clear()
    with pytest.raises(KeyError):
        gone.getNumberOfMatches()


def test_map_identification_data_view_edits_the_map():
    run, query, first, _, _ = make_run()
    features = oms.FeatureMap()
    view = features.identification_data_view().addRun(run)
    added = view.addMatch(query, run.getMatch(first), [0.2])
    stored = features.getIdentificationData().getRun("search")
    assert stored.getNumberOfMatches() == 3
    assert stored.getMatch(added).representation == "PEPTIDE"
    consensus = oms.ConsensusMap()
    consensus.identification_data_view().addRun("empty")
    assert [r.getIdentifier() for r in consensus.getIdentificationData().getRuns()] == ["empty"]


def test_value_types_accept_keyword_arguments_for_their_fields():
    parent = ID.QualifiedAccession(database="uniprot.fasta", accession="P02769")
    assert (parent.database, parent.accession) == ("uniprot.fasta", "P02769")
    assert parent == ID.QualifiedAccession(database="uniprot.fasta", accession="P02769")
    run = ID.Run("search")
    uniprot = run.addDatabase(ID.Database(path="uniprot.fasta", version="2026_04"))
    assert run.addDatabase(ID.Database(path="uniprot.fasta", version="2026_04")) == uniprot  # equal databases are one
    assert run.qualify(uniprot, "P02769") == parent and {uniprot: 1}[run.getDatabaseId(0)] == 1
    run.setDatabaseSequences([ID.DatabaseSequence(database=uniprot, accession="P02769", target_decoy=ID.TargetDecoy.TARGET)])
    assert run.getDatabaseSequences()[0].accession == "P02769"
    match = ID.MatchData(representation="LVNELTEFAK", charge=2,
                         sequence_evidence=[ID.SequenceEvidence(database=uniprot, accession="P02769", start=65, end=74,
                                                                before="K", after="T")])
    assert match.encoding == ID.Encoding.AA_SEQUENCE          # omitted fields keep their defaults
    assert match.sequence_evidence[0].accession == "P02769" and match.sequence_evidence[0].start == 65
    settings = ID.RunSettings(software="Comet", software_version="2024.01", metadata={"note": "test"})
    run.setSettings(settings)
    assert run.getSettings().software == "Comet" and run.getSettings().getMetaValue("note") == "test"
    search = oms.SearchParameters()
    search.db = "uniprot.fasta"
    with pytest.raises(Exception):  # the databases of a run are its own records
        run.setSettings(ID.RunSettings(software="Comet", search=search))
    compound = ID.MatchData(encoding=ID.Encoding.SMILES, representation="CCO", charge=1, formula="C2H6O",
                            identifiers=[ID.QualifiedAccession(database="HMDB", accession="HMDB0000108")],
                            adduct=oms.AdductInfo.parseAdductString("M+H;1+"), calculated_mz=47.0491)
    assert compound.adduct.getName() == "M+H;1+" and compound.identifiers[0].accession == "HMDB0000108"
    observation = ID.Observation(data_id="scan=1", rt=12.5, mz=None, metadata={"FWHM": 3.5, "label": "a"})
    assert observation.rt == 12.5 and observation.mz is None
    assert observation.getMetaValues() == {"FWHM": 3.5, "label": "a"}
    score = ID.ScoreDefinition(name="expect", higher_better=False, software="Comet", parameters={"tolerance": 10.0})
    assert (score.name, score.higher_better, score.software) == ("expect", False, "Comet")
    assert score.parameters.getMetaValue("tolerance") == 10.0
    # Feature links are values as well (the evidence of match names a database of this run).
    query = run.addIdentification(run.addSource(ID.SourceFile(path="a.mzML")), observation)
    match_id = run.addMatch(query, match)
    reference = ID.MatchReference(run_uuid=run.getUuid(), match=match_id)
    assert reference == ID.MatchReference(run_uuid=run.getUuid(), match=match_id)
    assert ID.MoleculeIdentity(encoding=ID.Encoding.SMILES, representation="CCO").representation == "CCO"
    assert File.Options(threads=2, replace_existing=True).threads == 2
    # Unknown names, wrong types and positional values are rejected.
    with pytest.raises(TypeError, match="unexpected keyword argument 'acession'"):
        ID.QualifiedAccession(database="db", acession="P1")
    with pytest.raises(TypeError, match="invalid type for field 'charge'"):
        ID.MatchData(charge="two")
    with pytest.raises(TypeError):
        ID.QualifiedAccession("db", "P1")
    with pytest.raises(TypeError, match="'metadata' expects a dict"):
        ID.Observation(metadata=[("FWHM", 3.5)])
    with pytest.raises(TypeError, match="invalid value for metadata 'FWHM'"):
        ID.Observation(metadata={"FWHM": object()})
    with pytest.raises(TypeError, match="unexpected keyword argument 'metadata'"):
        ID.QualifiedAccession(metadata={})                     # only records with metadata take it
    with pytest.raises(RuntimeError, match="unknown database"):
        ID.Run("other").addMatch(query, match)                 # evidence names a database of its own run


def test_values_compare_by_value_and_only_identities_hash():
    payload = ID.MatchData(representation="PEPTIDE", charge=2, metadata={"note": "x"})
    assert payload == copy.copy(payload) and not payload != ID.MatchData(payload)
    assert payload != ID.MatchData(representation="PEPTIDE", charge=3)
    assert payload != "PEPTIDE"
    # Mutable values compare by value and are unhashable; IDs, references and accessions are keys.
    for value in (payload, ID.ScoreDefinition(), ID.Observation(), ID.SourceFile(), ID.InferenceResult(), ID(),
                  ID.Run("search"), File.Options()):
        with pytest.raises(TypeError, match="unhashable"):
            hash(value)
    keys = {ID.QualifiedAccession(database="db", accession="P1"), ID.QualifiedAccession(database="db", accession="P1"),
            ID.QueryId(1), ID.QueryId(1), ID.MoleculeIdentity(encoding=ID.Encoding.SMILES, representation="CCO")}
    assert len(keys) == 3
    adducts = {oms.AdductInfo.parseAdductString("M+H;1+"), oms.AdductInfo.parseAdductString("M+H;1+")}
    assert len(adducts) == 1 and oms.AdductInfo is oms.AMSE_AdductInfo
    # Options without a C++ comparison compare field by field.
    assert File.Options() == File.Options() and File.Options(threads=2) != File.Options()
    assert File.ScanOptions(projection=File.Projection(all_scores=False)) != File.ScanOptions()

    # Records of a run compare by ID, payload, scores and candidates, not just their payload.
    run, query, first, second, score = make_run()
    assert run.getMatch(first) == run.getMatch(first)
    assert run.getMatch(first) != run.getMatch(second)
    before = run.getMatch(first)
    run.setScore(first, score, 0.5)
    assert before != run.getMatch(first) and before.getData() == run.getMatch(first).getData()
    assert run.getIdentification(query) == run.getIdentification(query)
    identification = run.getIdentification(query)
    run.setSelectedMatch(query, first)
    assert identification != run.getIdentification(query)
    assert run.getSources() == copy.copy(run).getSources()
    for record in (before, identification, run.getSources()[0]):
        with pytest.raises(TypeError, match="unhashable"):
            hash(record)
    data = ID()
    view = data.addRun("search")
    assert ID(data) == data and ID(data) != ID()
    # Views are handles: equal if they address the same run of the same dataset object.
    assert view == data.run_view("search") and len({view, data.run_view("search")}) == 1
    assert view != ID(data).run_view("search")


def test_values_have_readable_reprs():
    parent = ID.QualifiedAccession(database="uniprot.fasta", accession="P02769")
    assert repr(parent) == "IdentificationData.QualifiedAccession(database='uniprot.fasta', accession='P02769')"
    observation = ID.Observation(data_id="scan=1", rt=1.5, metadata={"label": "a"})
    assert repr(observation) == "IdentificationData.Observation(data_id='scan=1', rt=1.5, mz=None, metadata={'label': 'a'})"
    assert repr(ID.QueryId(3)) == "QueryId(3)" and repr(ID.MatchId(4)) == "MatchId(4)"
    assert repr(ID.Run("search").addScore(ID.ScoreDefinition(name="expect"))) == "<ScoreId 0>"  # bound to its run
    assert "adduct=AdductInfo(name='M+H;1+', formula='H1', charge=1, mol_multiplier=1)" in repr(
        ID.MatchData(encoding=ID.Encoding.SMILES, representation="CCO", adduct=oms.AdductInfo.parseAdductString("M+H;1+")))
    assert repr(File.Options(threads=2)).startswith("IdentificationDataFile.Options(batch_rows=")
    data = ID()
    view = data.addRun("search")
    query = view.addIdentification(view.addSource(ID.SourceFile()), observation)
    match = view.addMatch(query, ID.MatchData(representation="PEPTIDE"))
    assert repr(view) == f"RunView('search', uuid='{view.getUuid()}', kind=MoleculeKind.PEPTIDE, queries=1, matches=1)"
    assert repr(view.getMatch(match)).startswith("Match(id=1, scores=[], data=IdentificationData.MatchData(representation='PEPTIDE'")
    assert repr(view.getIdentification(query)).startswith("Identification(id=1, matches=1, selected=None, observation=")
    assert repr(data) == "IdentificationData(runs=['search'], inference_results=0)"
    uuid = view.getUuid()
    data.clear()
    assert repr(view) == f"RunView(uuid='{uuid}', removed)"


def test_nested_field_values_remain_owned_after_replacement():
    result = ID.InferenceResult()
    first_input = ID.InferenceInput()
    first_input.run_identifier = "first"
    result.inputs = [first_input]
    retained = result.inputs[0]
    result.inputs = []
    assert retained.run_identifier == "first"
    run, _, first, _, _ = make_run()
    evidence = run.getMatch(first).sequence_evidence[0]
    sequences = run.getMatch(first).sequence_evidence
    sequences.clear()
    assert evidence.accession == "P1"
    assert len(run.getMatch(first).sequence_evidence) == 1


def test_copy_preserves_ids_and_edits_are_independent():
    run, query, first, _, score = make_run()
    clone = copy.deepcopy(run)
    assert clone.getUuid() == run.getUuid()
    assert clone.getMatch(first).getId() == first
    clone.setScore(first, clone.getScoreId(score.value), 0.2)
    assert run.getScore(first, score) == 0.01
    assert clone.getScore(first, clone.getPrimaryScore()) == 0.2
    assert clone.getIdentification(query).data_id == "scan=1"


def test_callbacks_retain_safe_values_and_rollback_on_exceptions():
    run, _, first, _, score = make_run()
    retained = []

    def invalid_callback(match):
        retained.append(match)
        run.setScore(first, score, 0.5)
        return True

    with pytest.raises(Exception):
        run.filterMatches(invalid_callback)
    assert run.getScore(first, score) == 0.01
    assert run.getNumberOfMatches() == 2

    def throwing_callback(match):
        retained.append(match)
        if len(retained) > 2:
            raise ValueError("callback failed")
        return False

    with pytest.raises(ValueError, match="callback failed"):
        run.filterMatches(throwing_callback)
    assert run.getNumberOfMatches() == 2
    run.eraseMatches(lambda match: True)
    del run
    gc.collect()
    assert retained[0].representation == "PEPTIDE"
    assert retained[-1].getScores() == [0.1]


def test_transform_callbacks_use_owned_payloads_and_commit():
    run, _, first, _, _ = make_run()
    retained = []

    def transform(payload):
        retained.append(payload)
        payload.name = "updated"

    run.transformMatches(transform)
    assert run.getMatch(first).name == "updated"
    retained[0].name = "changed later"
    assert run.getMatch(first).name == "updated"
    run.eraseMatches(lambda match: True)
    assert retained[1].representation == "EDITPEP"


def test_filter_preserves_inference_provenance_or_explicitly_discards():
    data, query, first, second = make_data_with_inference()
    data.filterMatches(lambda match: match.getId() == first, ID.InferencePolicy.PRESERVE)
    run = data.getRun("search")
    assert run.findMatch(second) is None
    assert run.getIdentification(query).getSelectedMatch() is None
    result = data.getInferenceResults()[0]
    assert result.inputs[0].run_uuid == run.getUuid()
    assert result.inputs[0].selection == "all candidates"
    assert result.proteins.getHits()[0].getScore() == 0.99
    data.filterMatches(lambda match: True, ID.InferencePolicy.DISCARD)
    assert data.getInferenceResults() == []


@pytest.mark.parametrize("threads", [1, 4])
def test_native_round_trip_and_nonreused_ids(tmp_path, threads):
    data, query, first, second = make_data_with_inference()
    data.filterMatches(lambda match: match.getId() == first, ID.InferencePolicy.PRESERVE)
    path = str(tmp_path / "native")
    options = File.Options()
    options.threads = threads
    options.batch_rows = 1
    options.row_group_rows = 1
    File.store(path, data, options)
    loaded = File.load(path, options)
    run = loaded.getRun("search")
    assert run.getUuid() == data.getRun("search").getUuid()
    match = run.getMatch(first)
    assert match.representation == "PEPTIDE"
    assert match.metaValueExists("empty")
    assert match.getMetaValue("empty") is None
    assert match.getMetaValue("annotation") == ""
    assert match.getMetaValue("counts") == [1, 2, 3]
    assert run.getDatabase(match.sequence_evidence[0].database).path == "db.fasta"
    result = loaded.getInferenceResults()[0]
    assert result.inputs[0].run_uuid == run.getUuid()
    assert result.inputs[0].selection == "all candidates"
    added = run.addMatch(query, match, [0.05])
    assert added.value > second.value
    descriptor = File.inspect(path)[0]
    assert descriptor.match_count == 1
    assert descriptor.query_count == 1
    # matches.parquet names its score columns after the definitions
    assert descriptor.score_columns == ["score_posterior_error_probability"]
    assert File.loadRun(path, descriptor.uuid).getNumberOfMatches() == 1


@pytest.mark.parametrize("threads", [1, 4])
def test_native_projected_scan_keeps_callback_values_alive(tmp_path, threads):
    run, _, _, _, _ = make_run()
    data = ID()
    data.addRun(run)
    path = str(tmp_path / "native")
    options = File.Options()
    options.threads = threads
    options.batch_rows = 1
    options.row_group_rows = 1
    File.store(path, data, options)
    scan = File.ScanOptions()
    scan.buffering = options
    scan.runs = [run.getUuid()]
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
    assert batches[0][0].data.sequence_evidence == []
    assert queries[0].data.data_id == "scan=1"
    gc.collect()
    assert batches[0][0].match_id == 1


def search_result(rt_shift=0.0):
    """Five spectra, each with a target candidate and a lower-scoring decoy candidate."""
    data = ID()
    run = data.addRun("search")
    run.setPrimaryScore(run.addScore(ID.ScoreDefinition(name="hyperscore", higher_better=True)))
    source = run.addSource(ID.SourceFile(path="a.mzML"))
    for i, sequence in enumerate(["PEPTIDEA", "PEPTIDEC", "PEPTIDED", "PEPTIDEE", "PEPTIDEF"]):
        query = run.addIdentification(source, ID.Observation(data_id=f"scan={i}", rt=100.0 * (i + 1) + rt_shift, mz=500.0 + i))
        run.addMatch(query, ID.MatchData(representation=sequence, charge=2, target_decoy=ID.TargetDecoy.TARGET), [10.0 - i])
        run.addMatch(query, ID.MatchData(representation="DECOYPEPK", charge=2, target_decoy=ID.TargetDecoy.DECOY), [1.0 + i / 10])
    return data


def test_native_fdr_filters_alignment_and_mztab_export(tmp_path):
    data = search_result()
    hyperscore = data.getPrimaryScoreDefinition()
    qvalue = oms.FalseDiscoveryRate().applyToObservationMatches(data, hyperscore)
    assert qvalue.name == "PSM-level q-value" and not qvalue.higher_better
    assert [d.name for d in data.getScoreDefinitions()] == ["hyperscore", "PSM-level q-value"]
    run = data.getRuns()[0]
    q = run.bindScore(run.findScore(qvalue))
    values = [q(match) for block in run.getSources() for query in block.identifications for match in query.getMatches()]
    assert values == [0.0, None] * 5  # the best candidate per query; decoys get no value by default

    oms.IDFilter.keepBestMatchPerObservation(data, hyperscore)
    assert data.getRuns()[0].getNumberOfMatches() == 5
    oms.IDFilter.filterObservationMatchesByScore(data, hyperscore, 8.0)  # inclusive
    assert data.getRuns()[0].getNumberOfIdentifications() == 3
    with_decoys = search_result()
    oms.IDFilter.removeDecoys(with_decoys)
    assert with_decoys.getRuns()[0].getNumberOfMatches() == 5

    reference, shifted = search_result(), search_result(rt_shift=10.0)
    transformations = oms.MapAlignmentAlgorithmIdentification().align([reference, shifted], 0)
    assert len(transformations) == 2
    transformations[1].fitModel("linear")
    oms.MapAlignmentTransformer.transformRetentionTimes(shifted, transformations[1], True)
    assert shifted.getRuns()[0].getSources()[0].identifications[0].rt == pytest.approx(100.0)

    path = tmp_path / "search.mzTab"
    oms.MzTabFile().store(str(path), oms.IdentificationDataConverter.exportMzTab(reference))
    assert "PEPTIDEA" in path.read_text()


def test_converter_turns_feature_annotations_into_links_and_back():
    run = oms.ProteinIdentification()
    run.setIdentifier("search")
    peptide = oms.PeptideIdentification()
    peptide.setIdentifier("search")
    peptide.setScoreType("hyperscore")
    peptide.setRT(100.0)
    peptide.setMZ(500.0)
    peptide.setHits([oms.PeptideHit(7.0, 1, 2, oms.AASequence.fromString("PEPTIDE"))])
    feature = oms.Feature()
    feature.setPeptideIdentifications(oms.PeptideIdentificationList([peptide]))
    features = oms.FeatureMap()
    features.push_back(feature)
    features.setProteinIdentifications([run])

    oms.IdentificationDataConverter.importFeatureIDs(features)
    linked = features[0]
    assert len(linked.getPeptideIdentifications()) == 0 and features.getProteinIdentifications() == []
    (link,) = linked.getIDMatches()
    data = features.getIdentificationData()
    assert data.findRunByUuid(link.run_uuid).getMatch(link.match).representation == "PEPTIDE"

    oms.IdentificationDataConverter.exportFeatureIDs(features)
    restored = features[0]
    assert len(restored.getIDMatches()) == 0 and features.getIdentificationData().empty()
    assert restored.getPeptideIdentifications()[0].getHits()[0].getSequence().toString() == "PEPTIDE"


def test_native_filter_writes_reduced_dataset_with_explicit_inference_policy(tmp_path):
    data, query, first, second = make_data_with_inference()
    source = str(tmp_path / "source")
    destination = str(tmp_path / "reduced")
    File.store(source, data)
    File.filter(source, destination, lambda uuid, record: record.scores[0] < 0.05, ID.InferencePolicy.PRESERVE)
    result = File.load(destination)
    assert result.getRun("search").getNumberOfMatches() == 1
    assert result.getRun("search").getMatch(first).getId() == first
    assert result.getRun("search").getIdentification(query).getSelectedMatch() is None
    assert result.getInferenceResults()[0].inputs[0].run_uuid == data.getRun("search").getUuid()
    assert result.getInferenceResults()[0].inputs[0].selection == "all candidates"
    stripped = str(tmp_path / "stripped")
    File.filter(destination, stripped, lambda uuid, record: True, ID.InferencePolicy.DISCARD)
    assert File.load(stripped).getInferenceResults() == []
    with pytest.raises(Exception):
        File.store(source, data)


def test_failed_native_load_preserves_destination(tmp_path):
    run, _, first, _, _ = make_run()
    data = ID()
    data.addRun(run)
    path = str(tmp_path / "native")
    File.store(path, data)
    matches = list(Path(path).glob("matches.parquet"))
    assert len(matches) == 1
    matches[0].unlink()
    with pytest.raises(Exception):
        File.load(path, data)
    assert data.getRun("search").getMatch(first).representation == "PEPTIDE"
    assert data.getRun("search").getUuid() == run.getUuid()


def test_optional_adduct_and_empty_queries_round_trip(tmp_path):
    run = ID.Run("compounds", ID.MoleculeKind.COMPOUND)
    score = ID.ScoreDefinition()
    score.name = "similarity"
    run.setPrimaryScore(run.addScore(score))
    source = run.addSource(ID.SourceFile())
    query = run.addIdentification(source, ID.Observation())
    empty = run.addIdentification(source, ID.Observation())
    payload = ID.MatchData()
    payload.encoding = ID.Encoding.SMILES
    payload.representation = "CCO"
    payload.charge = 1
    payload.adduct = oms.AMSE_AdductInfo.parseAdductString("2M+H;1+")
    match_id = run.addMatch(query, payload, [0.9])
    data = ID()
    data.addRun(run)
    path = str(tmp_path / "native")
    File.store(path, data)
    loaded = File.load(path).getRun("compounds")
    assert loaded.getNumberOfIdentifications() == 2
    assert loaded.getIdentification(empty).getMatches() == []
    adduct = loaded.getMatch(match_id).adduct
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
    imported = oms.IdentificationDataAdapter.importLegacy([proteins], peptides)
    assert len(imported.queries) == 1
    values = imported.data
    assert values.getRuns()[0].getNumberOfMatches() == 1
    exported = oms.IdentificationDataAdapter.toLegacy(values)
    assert exported.losses == []
    assert exported.peptides[0].getHits()[0].getSequence().toString() == "PEPTIDE"
    assert exported.peptides[0].getHits()[0].getScore() == 0.03
    assert len(exported.queries) == 1


def test_inference_pools_runs_without_editing_candidates():
    data = ID()
    inputs = []
    for name in ("analysis-A", "analysis-B"):
        run, _, _, _, score = make_run(name)
        data.addRun(run)
        item = oms.IdentificationDataInference.Input()
        item.run_uuid = run.getUuid()
        item.score = score
        inputs.append(item)
    result = oms.IdentificationDataInference.infer(data, inputs, "pooled")
    assert len(result.inputs) == 2
    assert [entry.run_uuid for entry in result.inputs] == [entry.run_uuid for entry in inputs]
    assert len(result.proteins.getHits()) == 1
    assert [run.getNumberOfMatches() for run in data.getRuns()] == [2, 2]
    assert result.protein_score.scope == ID.ScoreScope.PROTEIN
    assert data.getInferenceResults() == []
    data.addInferenceResult(result)
    oms.IdentificationDataInference.retainProteins(result, [])
    assert result.proteins.getHits() == []
    assert len(data.getInferenceResults()[0].proteins.getHits()) == 1
    assert data.getRun("analysis-A").getNumberOfMatches() == 2


def test_file_handler_import_in_fresh_interpreter():
    # Module defaults must not depend on enum registration in a later module.
    subprocess.run(
        [sys.executable, "-c", "import pyopenms; assert pyopenms.FileHandler() is not None"],
        check=True,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("log_options", [{}, {"log": None}, {"log": oms.LogType.NONE}])
def test_file_handler_loads_and_stores_native_owning_values(tmp_path, log_options):
    run, _, first, _, _ = make_run()
    data = ID()
    data.addRun(run)
    path = str(tmp_path / "search.idparquet")
    handler = oms.FileHandler()
    handler.storeIdentificationData(path, data, **log_options)
    loaded = handler.loadIdentificationData(path, **log_options)
    assert loaded.getRun("search").getUuid() == run.getUuid()
    assert loaded.getRun("search").getMatch(first).representation == "PEPTIDE"


def test_psm_table_flattens_a_native_bundle(tmp_path):
    pa = pytest.importorskip("pyarrow")
    run, query, first, second, _ = make_run()
    other, _, _, _, _ = make_run("second search")
    data = ID()
    data.addRun(run)
    data.addRun(other)
    path = str(tmp_path / "search.idparquet")
    oms.FileHandler().storeIdentificationData(path, data)
    table = File.psm_table(path)
    assert isinstance(table, pa.Table)
    assert table.num_rows == 4
    rows = table.to_pylist()
    assert [row["run_identifier"] for row in rows] == ["search", "search", "second search", "second search"]
    row = rows[0]
    assert row["run_uuid"] == run.getUuid()
    assert row["reference_file_name"] == "/measurements/sample.mzML"
    assert row["query_id"] == query.value and row["match_id"] == first.value
    assert row["spectrum_reference"] == "scan=1"
    assert row["rt"] == 12.5 and row["observed_mz"] == 400.2
    assert row["peptidoform"] == "PEPTIDE" and row["precursor_charge"] == 2
    assert row["protein_accessions"] == ["P1"]
    assert row["target_decoy"] == "" and row["is_decoy"] is False
    assert row["score_type"] == "posterior error probability" and row["higher_score_better"] is False
    assert row["score"] == 0.01 and row["score_posterior_error_probability"] == 0.01
    assert rows[1]["match_id"] == second.value and rows[1]["peptidoform"] == "EDITPEP"
    only = File.psm_table(path, runs=["second search"])
    assert only.num_rows == 2 and set(only["run_identifier"].to_pylist()) == {"second search"}
    with pytest.raises(KeyError):
        File.psm_table(path, runs=["missing"])
    legacy = tmp_path / "legacy.idparquet"
    legacy.mkdir()
    (legacy / "psms.parquet").write_bytes(b"PAR1")
    with pytest.raises(ValueError, match="OpenMS 3.6"):
        File.psm_table(str(legacy))


def test_dataset_ordered_score_contract():
    first, _, _, _, _ = make_run("A")
    second, _, _, _, _ = make_run("B")
    extra = ID.ScoreDefinition()
    extra.name = "supplementary"
    first.addScore(extra)
    second.addScore(extra)
    data = ID()
    data.addRun(first)
    data.addRun(second)
    definition = data.getPrimaryScoreDefinition()
    data.setPrimaryScore(definition)
    data.validate()
    assert [score.name for score in data.getScoreDefinitions()] == [definition.name, "supplementary"]
    with pytest.raises(Exception):
        data.setPrimaryScore(extra)  # Null supplementary values cannot become primary.
    assert data.getPrimaryScoreDefinition().name == definition.name
    reordered = ID.Run("reordered")
    reordered.addScore(extra)
    reordered.setPrimaryScore(reordered.addScore(definition))
    with pytest.raises(Exception):
        data.addRun(reordered)
    incompatible, _, _, _, _ = make_run("bad")
    other = ID.ScoreDefinition()
    other.name = "different score"
    score_id = incompatible.addScore(other)
    for block in incompatible.getSources():
        for query in block.identifications:
            for match in query.getMatches():
                incompatible.setScore(match.getId(), score_id, 1.0)
    incompatible.setPrimaryScore(score_id)
    with pytest.raises(Exception):
        data.addRun(incompatible)
    assert len(data.getRuns()) == 2
    with pytest.raises(Exception):
        data.setPrimaryScore(other)
    assert data.getPrimaryScoreDefinition().name == definition.name


def test_merge_preserves_stable_links_with_repeated_display_names():
    first, query, match, _, _ = make_run()
    second, _, _, _, _ = make_run()
    data, other = ID(), ID()
    data.addRun(first)
    other.addRun(second)
    data.merge(other)
    assert len(data.getRuns()) == 2
    assert data.findRunByUuid(second.getUuid()).getIdentifier() == "search#2"
    reference = ID.MatchReference()
    reference.run_uuid = first.getUuid()
    reference.match = match
    assert data.findRunByUuid(reference.run_uuid).getMatch(reference.match).representation == "PEPTIDE"
    assert first.getIdentificationForMatch(match).getId() == query
    assert oms.IdentificationDataAdapter.QueryReference is ID.QueryReference
    duplicate = copy.copy(data)
    data.merge(duplicate)
    assert data == duplicate




def test_identity_handles_are_hashable_value_keys():
    run, query, first, second, _ = make_run()
    # Equal handles must hash equally so they work as dict keys and set members.
    assert hash(ID.MatchId(first.value)) == hash(first)
    scores = {first: 0.01, second: 0.1}
    assert scores[ID.MatchId(first.value)] == 0.01
    reference = ID.MatchReference()
    reference.run_uuid = run.getUuid()
    reference.match = first
    same = ID.MatchReference()
    same.run_uuid = run.getUuid()
    same.match = ID.MatchId(first.value)
    assert reference == same and hash(reference) == hash(same)
    assert len({reference, same}) == 1
    query_reference = ID.QueryReference()
    query_reference.run_uuid = run.getUuid()
    query_reference.query = query
    assert query_reference in {query_reference}


def test_predicates_use_python_truthiness():
    run, _, first, _, _ = make_run()
    # 0/1 and None are accepted like any Python condition, not only exact booleans.
    assert run.filterMatches(lambda match: 1 if match.getId() == first else 0) == 1
    assert run.getNumberOfMatches() == 1
    assert run.eraseMatches(lambda match: None) == 0


def test_native_destination_is_replaced_only_on_request(tmp_path):
    run, _, _, _, _ = make_run()
    data = ID()
    data.addRun(run)
    path = str(tmp_path / "search.idparquet")
    File.store(path, data)
    with pytest.raises(Exception):
        File.store(path, data)
    options = File.Options()
    options.replace_existing = True
    File.store(path, data, options)
    assert File.load(path) == data
    # Tools write through FileHandler, which replaces a previous native bundle like any other output.
    handler = oms.FileHandler()
    handler.storeIdentificationData(path, data)
    assert handler.loadIdentificationData(path) == data
    # Anything else at the destination is never replaced.
    other = tmp_path / "notes.idparquet"
    other.mkdir()
    (other / "keep.txt").write_text("user data")
    with pytest.raises(Exception):
        handler.storeIdentificationData(str(other), data)
    assert (other / "keep.txt").read_text() == "user data"


def test_derived_scores_record_their_producer():
    run = oms.ProteinIdentification()
    run.setSearchEngine("Comet")
    run.setSearchEngineVersion("2024.01")
    assert run.getScoreSoftware("expect") == ("Comet", "2024.01")
    run.setScoreSoftware("Posterior Error Probability", "IDPosteriorErrorProbability", "3.7.0")
    assert run.getScoreSoftware("Posterior Error Probability") == ("IDPosteriorErrorProbability", "3.7.0")
    run.clearScoreSoftware()
    assert run.getScoreSoftware("Posterior Error Probability") == ("Comet", "2024.01")
