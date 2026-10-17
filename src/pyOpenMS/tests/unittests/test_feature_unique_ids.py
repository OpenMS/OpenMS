# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
"""Every feature of a map needs a valid unique id that no other feature of the map shares."""

import pytest

import pyopenms as oms

INVALID = "without a valid unique id"
DUPLICATE = "Unique id shared by several"


def _feature(unique_id=None):
    feature = oms.Feature()
    if unique_id is not None:
        feature.setUniqueId(unique_id)
    return feature


def _feature_map(*unique_ids):
    fmap = oms.FeatureMap()
    for unique_id in unique_ids:
        fmap.push_back(_feature(unique_id))
    return fmap


def test_check_feature_map():
    oms.FeatureUniqueIds.check(oms.FeatureMap())
    oms.FeatureUniqueIds.check(_feature_map(1, 2))

    with pytest.raises(Exception, match=INVALID):
        oms.FeatureUniqueIds.check(_feature_map(1, None))
    with pytest.raises(Exception, match=DUPLICATE):
        oms.FeatureUniqueIds.check(_feature_map(1, 2, 2))

    # subordinates count too
    parent = _feature(1)
    parent.setSubordinates([_feature(2)])
    fmap = oms.FeatureMap()
    fmap.push_back(parent)
    fmap.push_back(_feature(3))
    oms.FeatureUniqueIds.check(fmap)

    parent.setSubordinates([_feature(3)])
    fmap = oms.FeatureMap()
    fmap.push_back(parent)
    fmap.push_back(_feature(3))
    with pytest.raises(Exception, match=DUPLICATE):
        oms.FeatureUniqueIds.check(fmap)


def test_check_consensus_map():
    def consensus_map(*unique_ids):
        cmap = oms.ConsensusMap()
        for unique_id in unique_ids:
            cf = oms.ConsensusFeature()
            if unique_id is not None:
                cf.setUniqueId(unique_id)
            cmap.push_back(cf)
        return cmap

    oms.FeatureUniqueIds.check(oms.ConsensusMap())
    oms.FeatureUniqueIds.check(consensus_map(1, 2))
    with pytest.raises(Exception, match=INVALID):
        oms.FeatureUniqueIds.check(consensus_map(1, None))
    with pytest.raises(Exception, match=DUPLICATE):
        oms.FeatureUniqueIds.check(consensus_map(1, 1))


def test_seed_list_generator_assigns_ids():
    seeds = [(1.1, 111.111), (2.2, 222.222)]
    features = oms.SeedListGenerator().convertSeedList(seeds)
    assert features.size() == 2
    oms.FeatureUniqueIds.check(features)
