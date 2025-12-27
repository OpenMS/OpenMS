#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Targeted Proteomics and MRM tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms
import copy
import numpy as np

# Import shared test helper functions
from test_helpers import report


@report
def testReactionMonitoringTransition():
    """
    @tests: ReactionMonitoringTransition
     ReactionMonitoringTransition.__init__
     ReactionMonitoringTransition.__eq__
     ReactionMonitoringTransition.__ne__
     ReactionMonitoringTransition.getName
     ReactionMonitoringTransition.setName
     ReactionMonitoringTransition.getNativeID
     ReactionMonitoringTransition.setNativeID
     ReactionMonitoringTransition.getPeptideRef
     ReactionMonitoringTransition.setPeptideRef
     ReactionMonitoringTransition.getCompoundRef
     ReactionMonitoringTransition.setCompoundRef
     ReactionMonitoringTransition.getProductMZ
     ReactionMonitoringTransition.setProductMZ
     ReactionMonitoringTransition.getPrecursorMZ
     ReactionMonitoringTransition.setPrecursorMZ
     ReactionMonitoringTransition.getLibraryIntensity
     ReactionMonitoringTransition.setLibraryIntensity
     ReactionMonitoringTransition.isDetectingTransition
     ReactionMonitoringTransition.setDetectingTransition
     ReactionMonitoringTransition.isIdentifyingTransition
     ReactionMonitoringTransition.setIdentifyingTransition
     ReactionMonitoringTransition.isQuantifyingTransition
     ReactionMonitoringTransition.setQuantifyingTransition
     ReactionMonitoringTransition.hasPrecursorCVTerms
     ReactionMonitoringTransition.setPrecursorCVTermList
     ReactionMonitoringTransition.addPrecursorCVTerm
     ReactionMonitoringTransition.getPrecursorCVTermList
     ReactionMonitoringTransition.addProductCVTerm
     ReactionMonitoringTransition.getIntermediateProducts
     ReactionMonitoringTransition.getProductChargeState
    """
    trans = pyopenms.ReactionMonitoringTransition()
    trans_ = copy.copy(trans)
    assert trans_ == trans
    trans_ = copy.deepcopy(trans)
    assert trans_ == trans
    trans_ = pyopenms.ReactionMonitoringTransition(trans)
    assert trans_ == trans

    trans.setName("name")
    assert trans.getName() == "name"
    trans.setNativeID("nativeID")
    assert trans.getNativeID() == "nativeID"
    trans.setPeptideRef("peptideRef")
    assert trans.getPeptideRef() == "peptideRef"
    trans.setCompoundRef("compoundRef")
    assert trans.getCompoundRef() == "compoundRef"

    assert isinstance(trans.getLibraryIntensity(), float)
    trans.setLibraryIntensity(17.0)
    assert abs(trans.getLibraryIntensity() - 17.0) < 1e-5

    assert isinstance(trans.getPrecursorMZ(), float)
    trans.setPrecursorMZ(5.0)
    assert trans.getPrecursorMZ() == 5.0

    assert isinstance(trans.getProductMZ(), float)
    trans.setProductMZ(4.0)
    assert trans.getProductMZ() == 4.0

    assert isinstance(trans.isDetectingTransition(), bool)
    trans.setDetectingTransition(False)
    assert not trans.isDetectingTransition()

    assert isinstance(trans.isIdentifyingTransition(), bool)
    trans.setIdentifyingTransition(True)
    assert trans.isIdentifyingTransition()

    assert isinstance(trans.isQuantifyingTransition(), bool)
    trans.setQuantifyingTransition(True)
    assert trans.isQuantifyingTransition()

    assert trans == trans
    assert not trans != trans

    # CV Terms
    assert not trans.hasPrecursorCVTerms()
    cvterm = pyopenms.CVTerm()
    trans.addPrecursorCVTerm(cvterm)
    assert trans.hasPrecursorCVTerms()

    trans.addProductCVTerm(cvterm)

    # Products - test basic interface without specific types
    assert len(trans.getIntermediateProducts()) == 0
    
    # Product charge state
    assert trans.getProductChargeState() == 0


@report
def testTargetedExperiment():
    """
    @tests: TargetedExperiment
     TargetedExperiment.__init__
     TargetedExperiment.clear
     TargetedExperiment.getCVs
     TargetedExperiment.setCVs
     TargetedExperiment.getTargetCVTerms
     TargetedExperiment.setTargetCVTerms
     TargetedExperiment.getPeptides
     TargetedExperiment.setPeptides
     TargetedExperiment.getProteins
     TargetedExperiment.setProteins
     TargetedExperiment.getTransitions
     TargetedExperiment.setTransitions
     TargetedExperiment.__eq__
     TargetedExperiment.__ne__
     """
    m = pyopenms.TargetedExperiment()
    m_ = copy.copy(m)
    assert m_ == m
    m_ = copy.deepcopy(m)
    assert m_ == m
    m_ = pyopenms.TargetedExperiment(m)
    assert m_ == m

    m.clear(True)
    m.setCVs(m.getCVs())

    targeted = m

    targeted.setCVs(targeted.getCVs())
    targeted.setTargetCVTerms(targeted.getTargetCVTerms())
    targeted.setPeptides(targeted.getPeptides())
    targeted.setProteins(targeted.getProteins())
    targeted.setTransitions(targeted.getTransitions())

    assert m == m
    assert not m != m


@report
def testTargetedExperimentHelper():
    """
    @tests: TargetedExperimentHelper
     TargetedExperimentHelper.RetentionTime
     RetentionTime.RTUnit
     RetentionTime.RTType
     RetentionTime.__init__
     RetentionTime.isRTset
     RetentionTime.setRT
    """
    rtu = pyopenms.RetentionTime.RTUnit()
    rtu = pyopenms.RetentionTime.RTUnit.SECOND
    rtu = pyopenms.RetentionTime.RTUnit.MINUTE
    rtt = pyopenms.RetentionTime.RTType()
    rtt = pyopenms.RetentionTime.RTType.LOCAL
    rtt = pyopenms.RetentionTime.RTType.NORMALIZED
    rtt = pyopenms.RetentionTime.RTType.IRT

    rt = pyopenms.RetentionTime()
    assert rt.software_ref is not None
    assert not rt.isRTset()
    rt.setRT(5.0)
    assert rt.isRTset()


@report
def testMRMTransitionGroup():
    """
    @tests: MRMTransitionGroup
     MRMTransitionGroup.__init__
     MRMTransitionGroup.setTransitionGroupID
     MRMTransitionGroup.getTransitionGroupID
     MRMTransitionGroup.getTransitions
     MRMTransitionGroup.addTransition
     MRMTransitionGroup.addChromatogram
     MRMTransitionGroup.addFeature
     MRMTransitionGroup.get_chromatogram_df
     MRMTransitionGroup.get_feature_df
    """
    mrmgroup = pyopenms.MRMTransitionGroupCP()
    assert mrmgroup is not None

    mrmgroup.setTransitionGroupID("this_id")
    assert mrmgroup.getTransitionGroupID() == "this_id"

    assert len(mrmgroup.getTransitions()) == 0
    mrmgroup.addTransition(pyopenms.ReactionMonitoringTransition(), "tr1")
    assert len(mrmgroup.getTransitions()) == 1

    # add data for testing df output
    ## test chromatogram df
    rt, intensity = [[1.0], [5]]
    chrom = pyopenms.MSChromatogram()
    chrom.set_peaks([rt, intensity])
    chrom.setNativeID("tr1")
    mrmgroup.addChromatogram(chrom, 'tr1')

    df = mrmgroup.get_chromatogram_df()
    # Default columns: rt, intensity, precursor_mz, precursor_charge, product_mz, native_id
    # chromatogram_type and comment are NOT included by default
    assert df.shape == (1, 6)
    assert df.loc[0, 'rt'] == 1.0
    assert df.loc[0, 'intensity'] == 5
    assert df.loc[0, 'native_id'] == 'tr1'

    # Test non-default columns (chromatogram_type, comment) via explicit selection
    df_all = mrmgroup.get_chromatogram_df(columns=['rt', 'intensity', 'chromatogram_type', 'comment'])
    assert df_all.shape == (1, 4)
    assert df_all.loc[0, 'chromatogram_type'] == 'MASS_CHROMATOGRAM'

    ## feature 1
    f1 = pyopenms.MRMFeature()
    f1.setRT(1.0)
    f1.setMetaValue(b'leftWidth', 0.5)
    f1.setMetaValue(b'rightWidth', 1.5)
    f1.setMetaValue(b'peak_apices_sum', 10.0)
    f1.setOverallQuality(0.5)
    f1.setUniqueId(1)
    f1.setIntensity(20.0)
    mrmgroup.addFeature(f1)

    # feature 2
    f2 = pyopenms.MRMFeature()
    f2.setRT(2.0)
    f2.setMetaValue(b'peak_apices_sum', 20.0)
    f2.setOverallQuality(1.0)
    f2.setUniqueId(2)
    f2.setIntensity(40.0)
    mrmgroup.addFeature(f2)

    df = mrmgroup.get_feature_df(meta_values=[b'leftWidth', b'rightWidth', b'peak_apices_sum'])
    assert df.shape == (2, 6)
    assert df.loc[1, 'leftWidth'] == 0.5
    assert df.loc[1, 'rightWidth'] == 1.5
    assert df.loc[1, 'peak_apices_sum'] == 10.0
    assert df.loc[1, 'intensity'] == 20.0
    assert df.loc[1, 'quality'] == 0.5
    assert df.loc[1, 'rt'] == 1.0

    assert np.isnan(df.loc[2, 'leftWidth'])
    assert np.isnan(df.loc[2, 'rightWidth'])
    assert df.loc[2, 'peak_apices_sum'] == 20.0
    assert df.loc[2, 'intensity'] == 40.0
    assert df.loc[2, 'quality'] == 1.0
    assert df.loc[2, 'rt'] == 2.0

    # If get "all" meta values should get the same result
    df = mrmgroup.get_feature_df(meta_values='all')
    assert df.shape == (2, 6)
    assert df.loc[1, 'leftWidth'] == 0.5
    assert df.loc[1, 'rightWidth'] == 1.5
    assert df.loc[1, 'peak_apices_sum'] == 10.0
    assert df.loc[1, 'intensity'] == 20.0
    assert df.loc[1, 'quality'] == 0.5
    assert df.loc[1, 'rt'] == 1.0

    assert np.isnan(df.loc[2, 'leftWidth'])
    assert np.isnan(df.loc[2, 'rightWidth'])
    assert df.loc[2, 'peak_apices_sum'] == 20.0
    assert df.loc[2, 'intensity'] == 40.0
    assert df.loc[2, 'quality'] == 1.0
    assert df.loc[2, 'rt'] == 2.0


@report
def testMRMMapping():
    """
    @tests: MRMMapping
     MRMMapping.__init__
     MRMMapping.mapExperiment
    """
    p = pyopenms.MRMMapping()
    assert p.mapExperiment is not None
    e = pyopenms.MSExperiment()
    c = pyopenms.MSChromatogram()
    e.addChromatogram(c)
    assert e.getNrChromatograms() == 1

    o = pyopenms.MSExperiment()
    t = pyopenms.TargetedExperiment()


@report
def testMRMDecoy():
    """
    @tests: MRMDecoy
      MRMDecoy.__init__
      MRMDecoy.generateDecoys
     """
    mrmdecoy = pyopenms.MRMDecoy()
    assert mrmdecoy is not None

    assert pyopenms.MRMDecoy().generateDecoys is not None


@report
def testMRMAssay():
    """
    @tests: MRMAssay
     MRMAssay.__init__
    """
    e = pyopenms.MRMAssay()
    assert e


@report
def testMRMIonSeries():
    """
    @tests: MRMIonSeries
     MRMIonSeries.__init__
    """
    e = pyopenms.MRMIonSeries()
    assert e


@report
def testConfidenceScoring():
    """
    @tests: ConfidenceScoring
      ConfidenceScoring.__init__
     """
    scoring = pyopenms.ConfidenceScoring()
    assert scoring is not None


if __name__ == "__main__":
    # Run all tests in this module
    import sys
    module = sys.modules[__name__]
    
    test_functions = [getattr(module, name) for name in dir(module) 
                     if name.startswith('test') and callable(getattr(module, name))]
    
    print(f"Running {len(test_functions)} targeted/MRM tests...")
    for test_func in test_functions:
        try:
            test_func()
            print(f"✅ {test_func.__name__}")
        except Exception as e:
            print(f"❌ {test_func.__name__}: {e}")
