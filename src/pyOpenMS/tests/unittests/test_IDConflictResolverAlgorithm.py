import pyopenms as oms


def _feature(uid, sequence=None, score=0.1):
    """A feature of one run, identified as `sequence` (unidentified if None)."""
    f = oms.Feature()
    f.setUniqueId(uid)
    if sequence is not None:
        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(sequence))
        hit.setScore(score)
        pep = oms.PeptideIdentification()
        pep.setHigherScoreBetter(False)
        pep.setHits([hit])
        ids = oms.PeptideIdentificationList()
        ids.append(pep)
        f.setPeptideIdentifications(ids)
    return f


def test_annotateTransfers():
    # run 0 identified, run 1 linked without an identification -> transferred to run 1
    transferred = oms.ConsensusFeature()
    transferred.insert(0, _feature(1, "PEPTIDEK"))
    transferred.insert(1, _feature(2))
    # identified in both runs -> no transfer
    identified = oms.ConsensusFeature()
    identified.insert(0, _feature(3, "PEPTIDER"))
    identified.insert(1, _feature(4, "PEPTIDER", 0.2))

    cmap = oms.ConsensusMap()
    cmap.push_back(transferred)
    cmap.push_back(identified)

    oms.IDConflictResolverAlgorithm.resolve(cmap, True)
    oms.IDConflictResolverAlgorithm.annotateTransfers(cmap)

    assert list(cmap[0].getMetaValue("mbr_transfer_map_indices")) == [1]
    assert not cmap[1].metaValueExists("mbr_transfer_map_indices")
