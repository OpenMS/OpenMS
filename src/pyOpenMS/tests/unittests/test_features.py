#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Feature detection and handling tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms
import copy
import numpy as np
import os
import pandas as pd

# Import shared test helper functions
from test_helpers import (
    report, 
    _testMetaInfoInterface,
    _testUniqueIdInterface,
    _testParam,
    _testStrOutput
)

# Helper function for testFeatureMap
def testDataProcessing(dp):
    """Helper function to test DataProcessing objects"""
    assert isinstance(dp, pyopenms.DataProcessing)

@report
def testFeatureFileOptions():
    """
    @tests: FeatureFileOptions
     FeatureFileOptions.__init__
     FeatureFileOptions.getLoadConvexHull
     FeatureFileOptions.getLoadSubordinates
     FeatureFileOptions.getMetadataOnly
     FeatureFileOptions.getSizeOnly
     FeatureFileOptions.setLoadConvexHull
     FeatureFileOptions.setLoadSubordinates
     FeatureFileOptions.setMetadataOnly
     FeatureFileOptions.setSizeOnly
    """

    fo = pyopenms.FeatureFileOptions()
    fo.getLoadConvexHull()
    fo.getLoadSubordinates()
    fo.getSizeOnly()
    assert fo.setLoadConvexHull is not None
    assert fo.setLoadSubordinates is not None
    assert fo.setMetadataOnly is not None
    assert fo.setSizeOnly is not None

@report
def testFeatureFinderAlgorithmPicked():
    """
    @tests: FeatureFinderAlgorithmPicked
     FeatureFinderAlgorithmPicked.__init__
     FeatureFinderAlgorithmPicked.getDefaults
     FeatureFinderAlgorithmPicked.getName
     FeatureFinderAlgorithmPicked.getParameters
     FeatureFinderAlgorithmPicked.getProductName
     FeatureFinderAlgorithmPicked.setName
     FeatureFinderAlgorithmPicked.setParameters
    """
    ff = pyopenms.FeatureFinderAlgorithmPicked()
    p = ff.getDefaults()
    _testParam(p)

    _testParam(ff.getParameters())

    assert ff.getName() == "FeatureFinderAlgorithmPicked"

    ff.setParameters(pyopenms.Param())

    ff.setName("test")
    assert ff.getName() == "test"


@report
def testFeatureDeconvolution():
    """
    @tests: FeatureDeconvolution
     FeatureDeconvolution.__init__
    """
    ff = pyopenms.FeatureDeconvolution()
    p = ff.getDefaults()
    _testParam(p)

    assert pyopenms.FeatureDeconvolution().compute is not None

@report
def testSeedListGenerator():
    """
    @tests: SeedListGenerator
     SeedListGenerator.__init__
    """
    ff = pyopenms.SeedListGenerator()

    # TODO
    # assert pyopenms.SeedListGenerator().generateSeedList is not None

    # TODO 
    # assert pyopenms.SeedListGenerator().compute is not None

@report
def testFeatureGrouping():
    """
    @tests: FeatureGroupingAlgorithmQT
     FeatureGroupingAlgorithmQT.getDefaults
     FeatureGroupingAlgorithmQT.getName
     FeatureGroupingAlgorithmQT.getParameters
     FeatureGroupingAlgorithmQT.setName
     FeatureGroupingAlgorithmQT.setParameters
     FeatureGroupingAlgorithmQT.transferSubelements
     FeatureGroupingAlgorithmQT.__init__
     FeatureGroupingAlgorithmQT.getDefaults
     FeatureGroupingAlgorithmQT.getName
     FeatureGroupingAlgorithmQT.getParameters
     FeatureGroupingAlgorithmQT.group
     FeatureGroupingAlgorithmQT.setName
     FeatureGroupingAlgorithmQT.setParameters
     FeatureGroupingAlgorithmQT.transferSubelements
    """

    assert pyopenms.FeatureGroupingAlgorithmQT.getDefaults is not None
    assert pyopenms.FeatureGroupingAlgorithmQT.getName is not None
    assert pyopenms.FeatureGroupingAlgorithmQT.getParameters is not None
    assert pyopenms.FeatureGroupingAlgorithmQT.setName is not None
    assert pyopenms.FeatureGroupingAlgorithmQT.setParameters is not None
    assert pyopenms.FeatureGroupingAlgorithmQT.transferSubelements is not None

    qt = pyopenms.FeatureGroupingAlgorithmQT()
    qt.getDefaults()
    qt.getParameters()
    qt.getName()
    assert qt.group is not None
    assert qt.setName is not None
    assert qt.setParameters is not None
    assert qt.transferSubelements is not None

@report
def testFeatureMap():
    """
    @tests: FeatureMap
     FeatureMap.__init__
     FeatureMap.__add__
     FeatureMap.__iadd__
     FeatureMap.__radd__
     FeatureMap.__getitem__
     FeatureMap.__iter__
     FeatureMap.__len__
     FeatureMap.append
     FeatureMap.clear
     FeatureMap.clearUniqueId
     FeatureMap.ensureUniqueId
     FeatureMap.extend
     FeatureMap.getDataProcessing
     FeatureMap.getProteinIdentifications
     FeatureMap.getUnassignedPeptideIdentifications
     FeatureMap.getUniqueId
     FeatureMap.setUniqueId
     FeatureMap.hasInvalidUniqueId
     FeatureMap.hasValidUniqueId
     FeatureMap.push_back
     FeatureMap.setDataProcessing
     FeatureMap.setProteinIdentifications
     FeatureMap.setUnassignedPeptideIdentifications
     FeatureMap.setUniqueIds
     FeatureMap.size
     FeatureMap.sortByIntensity
     FeatureMap.sortByMZ
     FeatureMap.sortByOverallQuality
     FeatureMap.sortByPosition
     FeatureMap.sortByRT
     FeatureMap.swap
     FeatureMap.updateRanges
    """
    fm = pyopenms.FeatureMap()
    fm_ = copy.copy(fm)
    assert fm_ == fm
    fm_ = copy.deepcopy(fm)
    assert fm_ == fm
    fm_ = pyopenms.FeatureMap(fm)
    assert fm_ == fm

    _testUniqueIdInterface(fm)
    fm.clear()
    fm.clearUniqueId()

    fm.getIdentifier()
    fm.getLoadedFileType()
    fm.getLoadedFilePath()

    f = pyopenms.Feature()
    fm.push_back(f)

    assert len(list(fm)) == 1


    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByIntensity()
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByIntensity(False)
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByPosition()
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByRT()
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByMZ()
    assert fm.size() == 1
    assert fm[0] == f

    fm.sortByOverallQuality()
    assert fm.size() == 1
    assert fm[0] == f

    fm2 = pyopenms.FeatureMap()

    fm.swap(fm2)
    assert fm2.size() == 1
    assert fm2[0] == f

    assert fm.size() == 0

    fm2.updateRanges()

    assert isinstance(fm2.getMinRT(), float)
    assert isinstance(fm2.getMinRT(), float)
    assert isinstance(fm2.getMaxMZ(), float)
    assert isinstance(fm2.getMaxMZ(), float)
    assert isinstance(fm2.getMinIntensity(), float)
    assert isinstance(fm2.getMaxIntensity(), float)

    assert fm2.getProteinIdentifications() == []
    fm2.setProteinIdentifications([])

    assert fm2.getUnassignedPeptideIdentifications().empty()
    fm2.setUnassignedPeptideIdentifications(pyopenms.PeptideIdentificationList())

    fm2.clear()
    assert fm2.size() == 0

    dp = pyopenms.DataProcessing()
    fm2.setDataProcessing([dp])
    assert fm2.getDataProcessing() == [dp]
    testDataProcessing(dp)

    fm2.setUniqueIds()

    fm += fm
    assert fm + fm2 != fm

    # Test __repr__ and __str__ methods
    fm_repr = pyopenms.FeatureMap()
    f1 = pyopenms.Feature()
    f1.setRT(100.0)
    f1.setMZ(500.0)
    f1.setIntensity(1000.0)
    f2 = pyopenms.Feature()
    f2.setRT(200.0)
    f2.setMZ(600.0)
    f2.setIntensity(2000.0)
    fm_repr.push_back(f1)
    fm_repr.push_back(f2)

    repr_str = repr(fm_repr)
    assert "FeatureMap(" in repr_str
    assert "num_features=" in repr_str

    # Test __len__, append, and extend methods
    fm_len = pyopenms.FeatureMap()
    assert len(fm_len) == 0
    assert len(fm_len) == fm_len.size()

    f_test1 = pyopenms.Feature()
    f_test1.setRT(100.0)
    f_test1.setMZ(500.0)

    f_test2 = pyopenms.Feature()
    f_test2.setRT(200.0)
    f_test2.setMZ(600.0)

    f_test3 = pyopenms.Feature()
    f_test3.setRT(300.0)
    f_test3.setMZ(700.0)

    # Test append (single item)
    fm_len.append(f_test1)
    assert len(fm_len) == 1
    assert len(fm_len) == fm_len.size()

    # Test extend with list
    fm_len.extend([f_test2, f_test3])
    assert len(fm_len) == 3
    assert len(fm_len) == fm_len.size()

    # Verify the features were added correctly
    assert fm_len[0].getRT() == 100.0
    assert fm_len[1].getRT() == 200.0
    assert fm_len[2].getRT() == 300.0

    # Test extend with another FeatureMap
    fm_source = pyopenms.FeatureMap()
    f_test4 = pyopenms.Feature()
    f_test4.setRT(400.0)
    f_test4.setMZ(800.0)
    fm_source.push_back(f_test4)

    fm_len.extend(fm_source)
    assert len(fm_len) == 4
    assert fm_len[3].getRT() == 400.0


@report
def testFeatureXMLFile():
    """
    @tests: FeatureXMLFile
     FeatureXMLFile.__init__
     FeatureXMLFile.load
     FeatureXMLFile.store
     FeatureXMLFile.getOptions
     FeatureXMLFile.setOptions
     FeatureXMLFile.loadSize

     FileHandler.__init__
     FileHandler.loadFeature
    """

    fm = pyopenms.FeatureMap()
    fm.setUniqueIds()

    f = pyopenms.Feature()
    f.setMZ(200)
    f.setCharge(1)
    f.setRT(10)
    f.setIntensity(10000)
    f.setOverallQuality(10)

    ch = pyopenms.ConvexHull2D()
    ch.setHullPoints(np.asarray([[8,199],[12,201]], dtype='f'))
    f.setConvexHulls([ch])

    f.setMetaValue(b'mv1', 1)
    f.setMetaValue(b'mv2', 2)

    f.setMetaValue('spectrum_native_id', 'spectrum=123')
    pep_id = pyopenms.PeptideIdentification()
    pep_id.insertHit(pyopenms.PeptideHit())
    peplist = pyopenms.PeptideIdentificationList()
    peplist.push_back(pep_id)
    f.setPeptideIdentifications(peplist)

    fm.push_back(f)

    f.setMetaValue('spectrum_native_id', 'spectrum=124')
    fm.push_back(f)

    assert fm.get_assigned_peptide_identifications().size() == 2
    assert fm.get_df(meta_values='all').shape == (2, 16)
    assert fm.get_df(meta_values='all', export_peptide_identifications=False).shape == (2, 12)

    assert pd.merge(fm.get_df(), pyopenms.peptide_identifications_to_df(fm.get_assigned_peptide_identifications()),
                on = ['feature_id', 'ID_native_id', 'ID_filename']).shape == (2,24)

    fm = pyopenms.FeatureMap()
    pyopenms.FeatureXMLFile().load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BSA1_F1_idmapped.featureXML'), fm)

    assert pd.merge(fm.get_df(), pyopenms.peptide_identifications_to_df(fm.get_assigned_peptide_identifications()),
                    on = ['feature_id', 'ID_native_id', 'ID_filename']).shape == (15,26)

    fh = pyopenms.FeatureXMLFile()
    fh.store("test.featureXML", fm)
    fh.load("test.featureXML", fm)

    fh = pyopenms.FileHandler()
    fh.loadFeatures("test.featureXML", fm)

@report
def testFileDescription():
    """
    @tests: ColumnHeader
     ColumnHeader.__init__
     ColumnHeader.filename
     ColumnHeader.label
     ColumnHeader.size
     ColumnHeader.unique_id
    """
    fd = pyopenms.ColumnHeader()
    _testStrOutput(fd.filename)
    _testStrOutput(fd.label)
    assert isinstance(fd.size, int)
    # assert isinstance(fd.unique_id, (long, int, bytes))

@report
def testElutionPeakDetection():
    """
    @tests: ElutionPeakDetection
     ElutionPeakDetection.__init__
     ElutionPeakDetection.detectPeaks
     ElutionPeakDetection.filterByPeakWidth
     ElutionPeakDetection.computeMassTraceNoise
     ElutionPeakDetection.computeMassTraceSNR
     ElutionPeakDetection.computeApexSNR
     ElutionPeakDetection.findLocalExtrema
     ElutionPeakDetection.smoothData
    """
    detection = pyopenms.ElutionPeakDetection()

    assert detection.detectPeaks is not None
    assert detection.filterByPeakWidth is not None
    assert detection.computeMassTraceNoise is not None
    assert detection.computeMassTraceSNR is not None
    assert detection.computeApexSNR is not None
    assert detection.findLocalExtrema is not None
    assert detection.smoothData is not None

    trace = pyopenms.Kernel_MassTrace()
    detection.smoothData(trace, 4)


if __name__ == "__main__":
    # Run all tests in this module
    import sys
    module = sys.modules[__name__]
    
    test_functions = [getattr(module, name) for name in dir(module) 
                     if name.startswith('test') and callable(getattr(module, name)) 
                     and name not in ['testDataProcessing']]  # Exclude helper function
    
    print(f"Running {len(test_functions)} feature tests...")
    for test_func in test_functions:
        try:
            test_func()
            print(f"✅ {test_func.__name__}")
        except Exception as e:
            print(f"❌ {test_func.__name__}: {e}")
