"""The filters of IDFilter for identification data and for the identifications of maps."""

import pytest
import pyopenms as oms


ID = oms.IdentificationData


def make_data():
    """A run "search" with q-values and proteins P1-P3: identification 1 (RT 10) has PEPTIDEB/P2 (0.05) and
    PEPTIDEA/P1 (0.01), identification 2 (RT 20) PEPTIDEC/P3 (0.02), the latter a decoy."""
    data = ID()
    run = data.addRun("search")
    score = run.addScore(ID.ScoreDefinition(name="q-value", higher_better=False))
    run.setPrimaryScore(score)
    database = run.addDatabase(ID.Database(path="proteins.fasta"))
    run.setDatabaseSequences([ID.DatabaseSequence(database=database, accession=accession) for accession in ("P1", "P2", "P3")])
    source = run.addSource(ID.SourceFile(path="sample.mzML"))
    queries = [run.addIdentification(source, ID.Observation(data_id="scan=%d" % i, rt=10.0 * i, mz=500.0)) for i in (1, 2)]

    def add(query, sequence, accession, value, decoy=False):
        evidence = ID.SequenceEvidence(database=database, accession=accession)
        match = ID.MatchData(representation=sequence, charge=2, sequence_evidence=[evidence],
                             target_decoy=ID.TargetDecoy.DECOY if decoy else ID.TargetDecoy.TARGET)
        return run.addMatch(query, match, [value])

    add(queries[0], "PEPTIDEB", "P2", 0.05)
    add(queries[0], "PEPTIDEA", "P1", 0.01)
    add(queries[1], "PEPTIDEC", "P3", 0.02, decoy=True)
    return data


def sequences(data):
    return [match.representation for source in data.getRuns()[0].getSources() for query in source.identifications for match in query.getMatches()]


def accessions(data):
    return [sequence.accession for sequence in data.getRuns()[0].getDatabaseSequences()]


def test_removeMatches():
    data = make_data()
    calls = []

    def remove(run, identification, match):
        calls.append(identification.data_id)
        # the predicate gets the live run, e.g. for the scores of the match
        return run.getScores(match)[0] > 0.03

    assert oms.IDFilter.removeMatches(data, remove) == 1
    assert sequences(data) == ["PEPTIDEA", "PEPTIDEC"]
    assert calls == ["scan=1", "scan=1", "scan=2"]


def test_removeMatches_failing_predicate_changes_nothing():
    data = make_data()

    def remove(run, identification, match):
        if match.representation == "PEPTIDEC":
            raise RuntimeError("failing predicate")
        return True

    with pytest.raises(RuntimeError, match="failing predicate"):
        oms.IDFilter.removeMatches(data, remove)
    assert sequences(data) == ["PEPTIDEB", "PEPTIDEA", "PEPTIDEC"]
    # the run cannot be edited from a predicate
    with pytest.raises(RuntimeError, match="Cannot modify a run"):
        oms.IDFilter.removeMatches(data, lambda run, identification, match: run.setPrimaryScore(None))
    assert sequences(data) == ["PEPTIDEB", "PEPTIDEA", "PEPTIDEC"]


def test_removeIdentifications():
    data = make_data()
    assert oms.IDFilter.removeIdentifications(data, lambda run, identification: identification.rt > 15) == 1
    assert sequences(data) == ["PEPTIDEB", "PEPTIDEA"]


def test_filters_of_identification_data():
    data = make_data()
    oms.IDFilter.filterPeptidesByRT(data, 0.0, 15.0)
    assert sequences(data) == ["PEPTIDEB", "PEPTIDEA"]

    data = make_data()
    oms.IDFilter.filterHitsByScore(data, 0.02)
    assert sequences(data) == ["PEPTIDEA", "PEPTIDEC"]

    # the best first, as PeptideIdentification.sort()
    data = make_data()
    oms.IDFilter.sortMatchesByScore(data)
    assert sequences(data) == ["PEPTIDEA", "PEPTIDEB", "PEPTIDEC"]
    oms.IDFilter.keepNBestPeptideHits(data, 1)
    assert sequences(data) == ["PEPTIDEA", "PEPTIDEC"]

    data = make_data()
    oms.IDFilter.removeDecoyHits(data)
    assert sequences(data) == ["PEPTIDEB", "PEPTIDEA"]
    oms.IDFilter.removeUnreferencedProteins(data)
    assert accessions(data) == ["P1", "P2"]

    data = make_data()
    oms.IDFilter.keepHitsMatchingProteins(data, {"P1", "P3"})
    assert sequences(data) == ["PEPTIDEA", "PEPTIDEC"]
    assert accessions(data) == ["P1", "P3"]

    data = make_data()
    oms.IDFilter.keepPeptidesWithMatchingSequences(data, oms.IDFilter.extractPeptideSequences(make_data()) - {"PEPTIDEB"})
    assert sequences(data) == ["PEPTIDEA", "PEPTIDEC"]
    oms.IDFilter.filterHitsByScore(data, 0.015)
    oms.IDFilter.removeEmptyIdentifications(data)
    assert data.getRuns()[0].getNumberOfIdentifications() == 1


def test_legacyMetaValue():
    data = make_data()
    run = data.run_view("search")
    match = run.getSources()[0].identifications[1].getMatches()[0]
    assert oms.IDFilter.legacyMetaValue(run, match, "target_decoy") == "decoy"
    assert oms.IDFilter.legacyMetaValue(run.copy(), match, "target_decoy") == "decoy"
    assert oms.IDFilter.legacyMetaValue(run, match, "missing") is None


def test_editProteinRuns():
    data = make_data()

    def edit(proteins):
        assert len(proteins) == 1
        hits = proteins[0].getHits()
        assert [hit.getAccession() for hit in hits] == ["P1", "P2", "P3"]
        proteins[0].setHits(hits[::-1][:2])

    oms.IDFilter.editProteinRuns(data, edit)
    assert accessions(data) == ["P3", "P2"]
    # proteins cannot be added
    with pytest.raises(RuntimeError, match="proteins can only be removed or reordered"):
        oms.IDFilter.editProteinRuns(data, lambda proteins: proteins[0].insertHit(oms.ProteinHit()))
    assert accessions(data) == ["P3", "P2"]


def test_filters_of_consensus_maps():
    cmap = oms.ConsensusMap()
    cmap.setIdentificationData(make_data())
    view = cmap.identification_data_view()
    run = view.getRuns()[0]
    queries = [query.getId() for query in run.getSources()[0].identifications]
    feature = oms.ConsensusFeature()
    feature.setUniqueId(1)
    feature.addIDQuery(ID.QueryReference(run_uuid=run.getUuid(), query=queries[0]))
    cmap.push_back(feature)

    # the predicate's run is a view of the map's identification data
    assert oms.IDFilter.removeIdentifications(cmap, lambda run, identification: run.getIdentifier() == "search" and identification.rt < 15) == 1
    # the link of the feature to the removed identification is gone
    assert len(cmap[0].getIDQueries()) == 0
    assert sequences(cmap.identification_data_view()) == ["PEPTIDEC"]
    oms.IDFilter.filterPeptidesByCharge(cmap, 3, 4)
    assert sequences(cmap.identification_data_view()) == []
