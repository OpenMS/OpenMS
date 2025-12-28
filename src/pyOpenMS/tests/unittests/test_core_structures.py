#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Core data structure tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms
import copy
import numpy as np

# Import shared test helper functions
from test_helpers import (
    report, 
    _testMetaInfoInterface,
    _testUniqueIdInterface,
    _testStrOutput
)

@report
def testPeak():
    """
    @tests: Peak1D
     Peak1D.__init__
     Peak1D.getIntensity
     Peak1D.getMZ
     Peak1D.setIntensity
     Peak1D.setMZ
     Peak1D.__eq__
     Peak1D.__ge__
     Peak1D.__gt__
     Peak1D.__le__
     Peak1D.__lt__
     Peak1D.__ne__
     Peak2D.__init__
     Peak2D.getIntensity
     Peak2D.getMZ
     Peak2D.getRT
     Peak2D.setIntensity
     Peak2D.setMZ
     Peak2D.setRT
     Peak2D.__eq__
     Peak2D.__ge__
     Peak2D.__gt__
     Peak2D.__le__
     Peak2D.__lt__
     Peak2D.__ne__
    """
    p1 = pyopenms.Peak1D()
    p1.setIntensity(12.0)
    assert p1.getIntensity() == 12.0
    p1.setMZ(13.0)
    assert p1.getMZ() == 13.0

    assert p1 == p1
    assert not p1 != p1

    p2 = pyopenms.Peak2D()
    assert p2 == p2
    assert not p2 != p2
    p2.setIntensity(22.0)
    assert p2.getIntensity() == 22.0
    p2.setMZ(23.0)
    assert p2.getMZ() == 23.0
    p2.setRT(45.0)
    assert p2.getRT() == 45.0

    # Test __repr__ and __str__ methods for Peak1D
    repr_str = repr(p1)
    assert "Peak1D(" in repr_str
    assert "mz=" in repr_str
    assert "intensity=" in repr_str
    str_str = str(p1)
    assert str_str == repr_str

    # Test __repr__ and __str__ methods for Peak2D
    repr_str = repr(p2)
    assert "Peak2D(" in repr_str
    assert "rt=" in repr_str
    assert "mz=" in repr_str
    assert "intensity=" in repr_str
    str_str = str(p2)
    assert str_str == repr_str


@report
def testChromatogramPeak():
    """
    @tests: ChromatogramPeak
     ChromatogramPeak.__init__
     ChromatogramPeak.__eq__
     ChromatogramPeak.__ge__
     ChromatogramPeak.__gt__
     ChromatogramPeak.__le__
     ChromatogramPeak.__lt__
     ChromatogramPeak.__ne__
     ChromatogramPeak.getIntensity
     ChromatogramPeak.getRT
     ChromatogramPeak.setIntensity
     ChromatogramPeak.setRT
    """
    p = pyopenms.ChromatogramPeak()
    assert p == p
    assert not p != p


    p.setIntensity(12.0)
    p.setRT(34.0)
    assert p.getIntensity() == 12.0
    assert p.getRT() == 34.0

    # Test __repr__ and __str__ methods
    repr_str = repr(p)
    assert "ChromatogramPeak(" in repr_str
    assert "rt=" in repr_str
    assert "intensity=" in repr_str
    str_str = str(p)
    assert str_str == repr_str

    # Test constructor with RT and intensity
    p2 = pyopenms.ChromatogramPeak(100.5, 5000.0)
    assert p2.getRT() == 100.5
    assert p2.getIntensity() == 5000.0

    # Test constructor with int RT (DPosition1 converter accepts int/float)
    p3 = pyopenms.ChromatogramPeak(200, 10000.0)
    assert p3.getRT() == 200.0
    assert p3.getIntensity() == 10000.0

    # Test getPosition returns float
    pos = p2.getPosition()
    assert isinstance(pos, float)
    assert pos == 100.5

    # Test setPosition with float
    p2.setPosition(300.0)
    assert p2.getRT() == 300.0
    assert p2.getPosition() == 300.0

    # Test setPosition with int (DPosition1 converter accepts int/float)
    p2.setPosition(400)
    assert p2.getRT() == 400.0


@report
def testMobilityPeak1DRepr():
    """
    @tests: MobilityPeak1D.__repr__
    """
    p = pyopenms.MobilityPeak1D()
    p.setMobility(1.234)
    p.setIntensity(10000.0)

    # Test __repr__ and __str__ methods
    repr_str = repr(p)
    assert "MobilityPeak1D(" in repr_str
    assert "mobility=" in repr_str
    assert "intensity=" in repr_str
    str_str = str(p)
    assert str_str == repr_str


@report
def testBaseFeature():
    """
    @tests: BaseFeature
     BaseFeature.__init__
     BaseFeature.clearUniqueId
     BaseFeature.ensureUniqueId
     BaseFeature.getCharge
     BaseFeature.getKeys
     BaseFeature.getMetaValue
     BaseFeature.getQuality
     BaseFeature.getUniqueId
     BaseFeature.getWidth
     BaseFeature.hasInvalidUniqueId
     BaseFeature.hasValidUniqueId
     BaseFeature.metaValueExists
     BaseFeature.removeMetaValue
     BaseFeature.setCharge
     BaseFeature.setMetaValue
     BaseFeature.setQuality
     BaseFeature.setWidth
     BaseFeature.clearMetaInfo
     BaseFeature.setUniqueId
    """
    bf = pyopenms.BaseFeature()
    _testMetaInfoInterface(bf)
    _testUniqueIdInterface(bf)
    bf.clearUniqueId()
    assert bf.ensureUniqueId()
    assert bf.getCharge() == 0
    assert isinstance(bf.getQuality(), float)
    assert isinstance(bf.getUniqueId(), (int, ))
    assert isinstance(bf.getWidth(), float)

    assert not bf.hasInvalidUniqueId()
    assert bf.hasValidUniqueId()


    _testMetaInfoInterface(bf)
    bf.setCharge(1)
    bf.setQuality(0.0)
    bf.setWidth(1.0)

@report
def testFeature():
    """
    @tests: Feature
     Feature.__init__
     Feature.clearUniqueId
     Feature.ensureUniqueId
     Feature.getCharge
     Feature.getIntensity
     Feature.getKeys
     Feature.getMZ
     Feature.getMetaValue
     Feature.getQuality
     Feature.getRT
     Feature.getUniqueId
     Feature.getWidth
     Feature.hasInvalidUniqueId
     Feature.hasValidUniqueId
     Feature.metaValueExists
     Feature.removeMetaValue
     Feature.setCharge
     Feature.setIntensity
     Feature.setMZ
     Feature.setMetaValue
     Feature.setQuality
     Feature.setRT
     Feature.setWidth
     Feature.__eq__
     Feature.__ge__
     Feature.__gt__
     Feature.__le__
     Feature.__lt__
     Feature.__ne__
     Feature.clearMetaInfo

     Feature.getConvexHulls
     Feature.getSubordinates
     Feature.setConvexHulls
     Feature.setSubordinates
     Feature.setUniqueId

     Feature.getPeptideIdentifications
     Feature.setPeptideIdentifications
    """
    f = pyopenms.Feature()
    _testMetaInfoInterface(f)
    _testUniqueIdInterface(f)

    f.setConvexHulls(f.getConvexHulls())
    f.setSubordinates(f.getSubordinates())
    f.setUniqueId(12345)

    assert f == f
    assert not f != f

    f.setCharge(-1)
    assert f.getCharge() == -1
    f.setIntensity(10.0)
    assert f.getIntensity() == 10.0
    f.setQuality(0, 20.0)
    assert f.getQuality(0) == 20.0
    f.setRT(30.0)
    assert f.getRT() == 30.0
    f.setWidth(40.0)
    assert f.getWidth() == 40.0

    p = f.getPeptideIdentifications()
    f.setPeptideIdentifications(p)

    # Test __repr__ and __str__ methods
    f_repr = pyopenms.Feature()
    f_repr.setRT(1234.5)
    f_repr.setMZ(445.678)
    f_repr.setIntensity(100000.0)
    f_repr.setCharge(2)
    f_repr.setOverallQuality(0.95)
    
    repr_str = repr(f_repr)
    assert "Feature(" in repr_str
    assert "rt=" in repr_str
    assert "mz=" in repr_str
    assert "intensity=" in repr_str
    assert "charge=" in repr_str
    assert "quality=" in repr_str
    
    # Test str method
    str_str = str(f_repr)
    assert str_str == repr_str


@report
def testConsensusFeature():
    """
    @tests: ConsensusFeature
     ConsensusFeature.__eq__
     ConsensusFeature.__ge__
     ConsensusFeature.__gt__
     ConsensusFeature.__le__
     ConsensusFeature.__lt__
     ConsensusFeature.__ne__
     ConsensusFeature.__init__
     ConsensusFeature.clearUniqueId
     ConsensusFeature.computeConsensus
     ConsensusFeature.computeDechargeConsensus
     ConsensusFeature.computeMonoisotopicConsensus
     ConsensusFeature.ensureUniqueId
     ConsensusFeature.getCharge
     ConsensusFeature.getKeys
     ConsensusFeature.getMetaValue
     ConsensusFeature.getQuality
     ConsensusFeature.getUniqueId
     ConsensusFeature.getWidth
     ConsensusFeature.hasInvalidUniqueId
     ConsensusFeature.hasValidUniqueId
     ConsensusFeature.insert
     ConsensusFeature.metaValueExists
     ConsensusFeature.removeMetaValue
     ConsensusFeature.setCharge
     ConsensusFeature.setMetaValue
     ConsensusFeature.setQuality
     ConsensusFeature.setWidth
     ConsensusFeature.clearMetaInfo
     ConsensusFeature.setUniqueId
     ConsensusFeature.size

     ConsensusFeature.getPeptideIdentifications
     ConsensusFeature.setPeptideIdentifications
    """


    f = pyopenms.ConsensusFeature()
    f_ = copy.copy(f)
    assert f_ == f
    f_ = copy.deepcopy(f)
    assert f_ == f
    f_ = pyopenms.ConsensusFeature(f)
    assert f_ == f

    _testUniqueIdInterface(f)
    _testMetaInfoInterface(f)

    f.setCharge(1)
    f.setQuality(2.0)
    f.setWidth(4.0)
    assert f.getCharge() == 1
    assert f.getQuality() == 2.0
    assert f.getWidth() == 4.0

    f.insert(0, pyopenms.Peak2D(), 1)
    f.insert(1, pyopenms.BaseFeature())
    f.insert(2, pyopenms.ConsensusFeature())

    f.computeConsensus()
    f.computeDechargeConsensus()
    f.computeMonoisotopicConsensus()

    assert f.size() >= 0

    p = f.getPeptideIdentifications()
    f.setPeptideIdentifications(p)

    # Test __repr__ and __str__ methods
    cf_repr = pyopenms.ConsensusFeature()
    cf_repr.setRT(1234.5)
    cf_repr.setMZ(445.678)
    cf_repr.setIntensity(100000.0)
    cf_repr.setCharge(2)
    cf_repr.insert(0, pyopenms.Peak2D(), 1)
    cf_repr.insert(1, pyopenms.BaseFeature())

    repr_str = repr(cf_repr)
    assert "ConsensusFeature(" in repr_str
    assert "rt=" in repr_str
    assert "mz=" in repr_str
    assert "intensity=" in repr_str
    assert "num_features=" in repr_str
    str_str = str(cf_repr)
    assert str_str == repr_str


@report
def testRichPeak():
    """
    @tests: RichPeak1D
     RichPeak1D.__init__
     RichPeak1D.getIntensity
     RichPeak1D.getKeys
     RichPeak1D.getMZ
     RichPeak1D.__eq__
     RichPeak1D.__ge__
     RichPeak1D.__gt__
     RichPeak1D.__le__
     RichPeak1D.__lt__
     RichPeak1D.__ne__
     RichPeak1D.getMetaValue
     RichPeak1D.clearMetaInfo
     RichPeak1D.isMetaEmpty
     RichPeak1D.metaValueExists
     RichPeak1D.removeMetaValue
     RichPeak1D.setIntensity
     RichPeak1D.setMZ
     RichPeak1D.setMetaValue
     RichPeak2D.__init__
     RichPeak2D.clearUniqueId
     RichPeak2D.clearMetaInfo
     RichPeak2D.isMetaEmpty
     RichPeak2D.ensureUniqueId
     RichPeak2D.getIntensity
     RichPeak2D.getKeys
     RichPeak2D.getMZ
     RichPeak2D.getMetaValue
     RichPeak2D.getRT
     RichPeak2D.getUniqueId
     RichPeak2D.hasInvalidUniqueId
     RichPeak2D.hasValidUniqueId
     RichPeak2D.metaValueExists
     RichPeak2D.removeMetaValue
     RichPeak2D.setIntensity
     RichPeak2D.setMZ
     RichPeak2D.setMetaValue
     RichPeak2D.setUniqueId
     RichPeak2D.setRT
     RichPeak2D.__eq__
     RichPeak2D.__ge__
     RichPeak2D.__gt__
     RichPeak2D.__le__
     RichPeak2D.__lt__
     RichPeak2D.__ne__
     """

    p2 = pyopenms.RichPeak2D()
    _testMetaInfoInterface(p2)
    _testUniqueIdInterface(p2)
    assert p2 == p2
    assert not p2 != p2
    p2.setMZ(22.0)
    p2.setIntensity(23.0)
    p2.setRT(43.0)
    assert p2.getMZ() == (22.0)
    assert p2.getIntensity() == (23.0)
    assert p2.getRT() == (43.0)


@report
def testDPosition():
    dp = pyopenms.DPosition1()
    dp = pyopenms.DPosition1(1.0)
    assert dp[0] == 1.0

    dp = pyopenms.DPosition2()
    dp = pyopenms.DPosition2(1.0, 2.0)

    assert dp[0] == 1.0
    assert dp[1] == 2.0

    # Test __repr__ method for DPosition1
    dp1_repr = pyopenms.DPosition1(123.456)
    repr_str = repr(dp1_repr)
    assert "DPosition1(" in repr_str
    assert "x=" in repr_str

    # Test __repr__ method for DPosition2
    dp2_repr = pyopenms.DPosition2(100.0, 200.0)
    repr_str = repr(dp2_repr)
    assert "DPosition2(" in repr_str
    assert "x=" in repr_str
    assert "y=" in repr_str

@report
def testConvexHull2D():
    """
    @tests: ConvexHull2D
     ConvexHull2D.__eq__
     ConvexHull2D.__ge__
     ConvexHull2D.__gt__
     ConvexHull2D.__init__
     ConvexHull2D.__le__
     ConvexHull2D.__lt__
     ConvexHull2D.__ne__
     ConvexHull2D.clear
     """
    ch = pyopenms.ConvexHull2D()
    ch.clear()
    assert ch == ch

    # Test __repr__ and __str__ methods
    ch_repr = pyopenms.ConvexHull2D()
    ch_repr.addPointXY(100.0, 400.0)
    ch_repr.addPointXY(110.0, 400.0)
    ch_repr.addPointXY(110.0, 410.0)
    ch_repr.addPointXY(100.0, 410.0)

    repr_str = repr(ch_repr)
    assert "ConvexHull2D(" in repr_str
    assert "num_points=" in repr_str


@report
def testMRMFeature():
    """
    @tests: MRMFeature
      MRMFeature.__init__
      MRMFeature.addScore
      MRMFeature.getScore
     """
    mrmfeature = pyopenms.MRMFeature()
    f = pyopenms.Feature()

    fs = mrmfeature.getFeatures()
    assert len(fs) == 0

    mrmfeature.addFeature(f, "myFeature")
    fs = mrmfeature.getFeatures()
    assert len(fs) == 1
    assert mrmfeature.getFeature("myFeature") is not None
    slist = []
    mrmfeature.getFeatureIDs(slist)
    assert len(slist) == 1

    mrmfeature.addPrecursorFeature(f, "myFeature_Pr0")
    slist = []
    mrmfeature.getPrecursorFeatureIDs(slist)
    assert len(slist) == 1
    assert mrmfeature.getPrecursorFeature("myFeature_Pr0") is not None

    s = mrmfeature.getScores()
    assert abs(s.yseries_score - 0.0) < 1e-4
    s.yseries_score = 4.0
    mrmfeature.setScores(s)
    s2 = mrmfeature.getScores()
    assert abs(s2.yseries_score - 4.0) < 1e-4

    # Test __repr__ and __str__ methods
    mrm_repr = pyopenms.MRMFeature()
    mrm_repr.setRT(1234.5)
    mrm_repr.setMZ(445.678)
    mrm_repr.setIntensity(100000.0)
    mrm_repr.setCharge(2)
    f1 = pyopenms.Feature()
    f2 = pyopenms.Feature()
    mrm_repr.addFeature(f1, "trans1")
    mrm_repr.addFeature(f2, "trans2")

    repr_str = repr(mrm_repr)
    assert "MRMFeature(" in repr_str
    assert "rt=" in repr_str
    assert "mz=" in repr_str
    assert "intensity=" in repr_str
    assert "num_transitions=" in repr_str
    str_str = str(mrm_repr)
    assert str_str == repr_str


@report
def testAnnotationState():
    """
    @tests: AnnotationState
     AnnotationState.__init__
    """
    state = pyopenms.AnnotationState()

    assert state.FEATURE_ID_NONE is not None
    assert state.FEATURE_ID_SINGLE is not None
    assert state.FEATURE_ID_MULTIPLE_SAME is not None
    assert state.FEATURE_ID_MULTIPLE_DIVERGENT is not None
    assert state.SIZE_OF_ANNOTATIONSTATE is not None

@report
def testStringDataArray():
    """
    @tests: StringDataArray
     """
    da = pyopenms.StringDataArray()
    assert da.size() == 0
    da.push_back("hello")
    da.push_back("world")
    assert da.size() == 2
    assert da[0] == b"hello"
    assert da[1] == b"world"
    da[1] = b"hello world"
    assert da[1] == b"hello world", da[1]
    da.clear()
    assert da.size() == 0
    da.push_back("hello")
    assert da.size() == 1
    da.resize(3)
    da[0] = b"hello"
    da[1] = b""
    da[2] = b"world"
    assert da.size() == 3

    # Test __repr__ method
    sda_repr = pyopenms.StringDataArray()
    sda_repr.setName("annotation")
    sda_repr.push_back("test1")
    sda_repr.push_back("test2")

    repr_str = repr(sda_repr)
    assert "StringDataArray(" in repr_str
    assert "name='annotation'" in repr_str
    assert "size=2" in repr_str

@report
def testIntegerDataArray():
    """
    @tests: IntegerDataArray
     """
    da = pyopenms.IntegerDataArray()
    assert da.size() == 0
    da.push_back(1)
    da.push_back(4)
    assert da.size() == 2
    assert da[0] == 1
    assert da[1] == 4
    da[1] = 7
    assert da[1] == 7
    da.clear()
    assert da.size() == 0
    da.push_back(1)
    assert da.size() == 1
    da.resize(3)
    da[0] = 1
    da[1] = 2
    da[2] = 3
    assert da.size() == 3

    q = da.get_data()
    q = np.append(q, 4).astype(np.intc)
    da.set_data(q)
    assert da.size() == 4

    # Test __repr__ method
    ida_repr = pyopenms.IntegerDataArray()
    ida_repr.setName("charge_state")
    ida_repr.push_back(1)
    ida_repr.push_back(2)
    ida_repr.push_back(3)

    repr_str = repr(ida_repr)
    assert "IntegerDataArray(" in repr_str
    assert "name='charge_state'" in repr_str
    assert "size=3" in repr_str

@report
def testFloatDataArray():
    """
    @tests: FloatDataArray
     """
    da = pyopenms.FloatDataArray()
    assert da.size() == 0
    da.push_back(1.0)
    da.push_back(4.0)
    assert da.size() == 2
    assert da[0] == 1.0
    assert da[1] == 4.0
    da[1] = 7.0
    assert da[1] == 7.0
    da.clear()
    assert da.size() == 0
    da.push_back(1.0)
    assert da.size() == 1
    da.resize(3)
    da[0] = 1.0
    da[1] = 2.0
    da[2] = 3.0
    assert da.size() == 3

    q = da.get_data()
    q = np.append(q, 4.0).astype(np.float32)
    da.set_data(q)
    assert da.size() == 4

    # Test __repr__ method
    fda_repr = pyopenms.FloatDataArray()
    fda_repr.setName("ion_mobility")
    fda_repr.push_back(1.5)
    fda_repr.push_back(2.5)

    repr_str = repr(fda_repr)
    assert "FloatDataArray(" in repr_str
    assert "name='ion_mobility'" in repr_str
    assert "size=2" in repr_str

@report
def testGaussFitResult():
    """
    @tests: GaussFitResult
     GaussFitResult.__init__
    """
    ins = pyopenms.GaussFitResult(0.0, 0.0, 0.0)
    ins.A = 5.0
    ins.x0 = 5.0
    ins.sigma = 5.0

@report
def testGaussFitter():
    """
    @tests: GaussFitter
     GaussFitter.__init__
    """
    ins = pyopenms.GaussFitter()

@report
def testKernelMassTrace():
    trace = pyopenms.Kernel_MassTrace()

    assert trace.getSize is not None
    assert trace.getLabel is not None
    assert trace.setLabel is not None

    assert trace.getCentroidMZ is not None
    assert trace.getCentroidRT is not None
    assert trace.getCentroidSD is not None
    assert trace.getFWHM is not None
    assert trace.getTraceLength is not None
    assert trace.getFWHMborders is not None
    assert trace.getSmoothedIntensities is not None
    assert trace.getAverageMS1CycleTime is not None

    assert trace.computeSmoothedPeakArea is not None
    assert trace.computePeakArea is not None
    assert trace.findMaxByIntPeak is not None
    assert trace.estimateFWHM is not None
    assert trace.computeFwhmArea is not None
    assert trace.computeFwhmAreaSmooth is not None
    # assert trace.computeFwhmAreaRobust is not None
    # assert trace.computeFwhmAreaSmoothRobust is not None
    assert trace.getIntensity is not None
    assert trace.getMaxIntensity is not None

    assert trace.getConvexhull is not None

    assert trace.setCentroidSD is not None
    assert trace.setSmoothedIntensities is not None
    assert trace.updateSmoothedMaxRT is not None
    assert trace.updateWeightedMeanRT is not None
    assert trace.updateSmoothedWeightedMeanRT is not None
    assert trace.updateMedianRT is not None
    assert trace.updateMedianMZ is not None
    assert trace.updateMeanMZ is not None
    assert trace.updateWeightedMeanMZ is not None
    assert trace.updateWeightedMZsd is not None

    s = trace.getSize()


if __name__ == "__main__":
    # Run all tests in this module
    import sys
    module = sys.modules[__name__]
    
    test_functions = [getattr(module, name) for name in dir(module) 
                     if name.startswith('test') and callable(getattr(module, name))]
    
    print(f"Running {len(test_functions)} core structure tests...")
    for test_func in test_functions:
        try:
            test_func()
            print(f"✅ {test_func.__name__}")
        except Exception as e:
            print(f"❌ {test_func.__name__}: {e}")
