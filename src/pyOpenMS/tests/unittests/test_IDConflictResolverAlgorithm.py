import pyopenms as oms


def _feature(uid):
    """A feature of one run."""
    f = oms.Feature()
    f.setUniqueId(uid)
    return f


def _identification(map_index, sequence, score=0.1):
    """A peptide identification of run `map_index` (its "map_index", as a feature linker annotates it)."""
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString(sequence))
    hit.setScore(score)
    pep = oms.PeptideIdentification()
    pep.setIdentifier("search")
    pep.setScoreType("Posterior Error Probability")
    pep.setHigherScoreBetter(False)
    pep.setHits([hit])
    pep.setMetaValue("map_index", map_index)
    return pep


def _identifications(*peps):
    ids = oms.PeptideIdentificationList()
    for pep in peps:
        ids.append(pep)
    return ids


def test_annotateTransfers():
    # run 0 identified, run 1 linked without an identification -> transferred to run 1
    transferred = oms.ConsensusFeature()
    transferred.insert(0, _feature(1))
    transferred.insert(1, _feature(2))
    # identified in both runs -> no transfer
    identified = oms.ConsensusFeature()
    identified.insert(0, _feature(3))
    identified.insert(1, _feature(4))

    cmap = oms.ConsensusMap()
    cmap.push_back(transferred)
    cmap.push_back(identified)
    run = oms.ProteinIdentification()
    run.setIdentifier("search")
    ids = oms.MapPeptideIdentifications()
    ids.proteins = [run]
    ids.setFeature([0], _identifications(_identification(0, "PEPTIDEK")))
    ids.setFeature([1], _identifications(_identification(0, "PEPTIDER"), _identification(1, "PEPTIDER", 0.2)))
    oms.IdentificationDataConverter.importMapIdentifications(cmap, ids)

    oms.IDConflictResolverAlgorithm.resolve(cmap, True)
    oms.IDConflictResolverAlgorithm.annotateTransfers(cmap)

    assert list(cmap[0].getMetaValue("mbr_transfer_map_indices")) == [1]
    assert not cmap[1].metaValueExists("mbr_transfer_map_indices")
