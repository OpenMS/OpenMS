#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## File I/O tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms
import os

# Import shared test helper functions
from test_helpers import report


@report
def testFileHandler():
  """
  @tests: FileHandler
   FileHandler.__init__
   FileHandler.getType
   FileHandler.loadExperiment
   FileHandler.storeExperiment
  """
  mse = pyopenms.MSExperiment()

  fh = pyopenms.FileHandler()
  fh.storeExperiment("test1.mzML", mse)
  fh.loadExperiment("test1.mzML", mse)
  fh.storeExperiment("test1.mzXML", mse)
  fh.loadExperiment("test1.mzXML", mse)
  fh.storeExperiment("test1.mzData", mse)
  fh.loadExperiment("test1.mzData", mse)


@report
def testCachedMzML():
  """
  @tests: CachedmzML
   CachedmzML.store
   CachedmzML.load
   CachedmzML.getMetaData
   CachedmzML.getNrChromatograms
   CachedmzML.getNrSpectra
  """
  mse = pyopenms.MSExperiment()
  s = pyopenms.MSSpectrum()
  mse.addSpectrum(s)

  # First load data and cache to disk
  pyopenms.CachedmzML.store("myCache.mzML", mse)

  # Now load data
  cfile = pyopenms.CachedmzML()
  pyopenms.CachedmzML.load("myCache.mzML", cfile)

  meta_data = cfile.getMetaData()
  assert cfile.getNrChromatograms() == 0
  assert cfile.getNrSpectra() == 1


@report
def testIndexedMzMLFile():
  """
  @tests: IndexedMzMLHandler
   IndexedMzMLHandler.__init__
   IndexedMzMLHandler.getNrChromatograms
   IndexedMzMLHandler.getNrSpectra
   IndexedMzMLHandler.getMSSpectrumById
   IndexedMzMLHandler.getSpectrumById
  """
  mse = pyopenms.MSExperiment()
  s = pyopenms.MSSpectrum()
  mse.addSpectrum(s)

  # First load data and cache to disk
  pyopenms.MzMLFile().store("tfile_idx.mzML", mse)

  # Now load data
  ih = pyopenms.IndexedMzMLHandler("tfile_idx.mzML")

  assert ih.getNrChromatograms() == 0
  assert ih.getNrSpectra() == 1

  s = ih.getMSSpectrumById(0)
  s2 = ih.getSpectrumById(0)


@report
def testIndexedMzMLDecoder():
  """
  @tests: IndexedMzMLDecoder
   IndexedMzMLDecoder.__init__
   IndexedMzMLDecoder.findIndexListOffset
  """
  decoder = pyopenms.IndexedMzMLDecoder()

  try:
    pos = decoder.findIndexListOffset("abcde", 100)
    raise Exception("Should raise an error")
  except RuntimeError:
    pass


@report
def testConsensusXMLFile():
  """
  @tests: ConsensusXMLFile
   ConsensusXMLFile.__init__
   ConsensusXMLFile.getOptions
   ConsensusXMLFile.load
   ConsensusXMLFile.store
  """
  f = pyopenms.ConsensusXMLFile()
  f.getOptions()
  assert f.load is not None
  assert f.store is not None


@report
def testIdXMLFile():
  """
  @tests: IdXMLFile
   IdXMLFile.__init__
   IdXMLFile.load
   IdXMLFile.store
  """
  assert pyopenms.IdXMLFile().load is not None
  assert pyopenms.IdXMLFile().store is not None


@report
def testXTandemXMLFile():
  """
  @tests: XTandemXMLFile
   XTandemXMLFile.__init__
   XTandemXMLFile.load
   XTandemXMLFile.setModificationDefinitionsSet
  """
  f = pyopenms.XTandemXMLFile()
  assert f.load is not None


@report
def testXTandemInfile():
  """
  @tests: XTandemInfile
   XTandemInfile.__init__
   XTandemInfile.setFragmentMassTolerance
   XTandemInfile.getFragmentMassTolerance
   XTandemInfile.setPrecursorMassTolerancePlus
   XTandemInfile.getPrecursorMassTolerancePlus
   XTandemInfile.setPrecursorMassToleranceMinus
   XTandemInfile.getPrecursorMassToleranceMinus
   XTandemInfile.setPrecursorErrorType
   XTandemInfile.getPrecursorErrorType
   XTandemInfile.setFragmentMassErrorUnit
   XTandemInfile.getFragmentMassErrorUnit
   XTandemInfile.setPrecursorMassErrorUnit
   XTandemInfile.getPrecursorMassErrorUnit
  """
  f = pyopenms.XTandemInfile()

  assert f.setFragmentMassTolerance is not None
  assert f.getFragmentMassTolerance is not None

  assert f.setPrecursorMassTolerancePlus is not None
  assert f.getPrecursorMassTolerancePlus is not None
  assert f.setPrecursorMassToleranceMinus is not None
  assert f.getPrecursorMassToleranceMinus is not None

  assert f.setPrecursorErrorType is not None
  assert f.getPrecursorErrorType is not None

  assert f.setFragmentMassErrorUnit is not None
  assert f.getFragmentMassErrorUnit is not None
  assert f.setPrecursorMassErrorUnit is not None
  assert f.getPrecursorMassErrorUnit is not None


@report
def testPepXMLFile():
  """
  @tests: PepXMLFile
   PepXMLFile.__init__
   PepXMLFile.load
   PepXMLFile.store
  """
  f = pyopenms.PepXMLFile()

  assert pyopenms.PepXMLFile().load is not None
  assert pyopenms.PepXMLFile().store is not None


@report
def testProtXMLFile():
  """
  @tests: ProtXMLFile
   ProtXMLFile.__init__
   ProtXMLFile.load
   ProtXMLFile.store
  """
  f = pyopenms.ProtXMLFile()

  assert pyopenms.ProtXMLFile().load is not None
  assert pyopenms.ProtXMLFile().store is not None


@report
def testMzIdentMLFile():
  """
  @tests: MzIdentMLFile
   MzIdentMLFile.__init__
   MzIdentMLFile.load
   MzIdentMLFile.store
   MzIdentMLFile.isSemanticallyValid
  """
  f = pyopenms.MzIdentMLFile()

  assert pyopenms.MzIdentMLFile().load is not None
  assert pyopenms.MzIdentMLFile().store is not None
  assert pyopenms.MzIdentMLFile().isSemanticallyValid is not None


@report
def testMzTabFile():
  """
  @tests: MzTabFile
   MzTabFile.__init__
  """
  f = pyopenms.MzTabFile()


@report
def testMzTab():
  """
  @tests: MzTab
   MzTab.__init__
  """
  # Placeholder - MzTab class may not be fully implemented yet
  pass


@report
def testDTA2DFile():
  """
  @tests: DTA2DFile
   DTA2DFile.__init__
   DTA2DFile.load
   DTA2DFile.store
  """
  f = pyopenms.DTA2DFile()

  assert pyopenms.DTA2DFile().load is not None
  assert pyopenms.DTA2DFile().store is not None


@report
def testDTAFile():
  """
  @tests: DTAFile
   DTAFile.__init__
   DTAFile.load
   DTAFile.store
  """
  f = pyopenms.DTAFile()

  assert pyopenms.DTAFile().load is not None
  assert pyopenms.DTAFile().store is not None


@report
def testEDTAFile():
  """
  @tests: EDTAFile
   EDTAFile.__init__
   EDTAFile.load
   EDTAFile.store
  """
  f = pyopenms.EDTAFile()

  assert pyopenms.EDTAFile().load is not None
  assert pyopenms.EDTAFile().store is not None


@report
def testKroenikFile():
  """
  @tests: KroenikFile
   KroenikFile.__init__
   KroenikFile.load
   KroenikFile.store
  """
  f = pyopenms.KroenikFile()

  assert pyopenms.KroenikFile().load is not None
  assert pyopenms.KroenikFile().store is not None


@report
def testMSPFile():
  """
  @tests: MSPFile
   MSPFile.__init__
  """
  f = pyopenms.MSPFile()


@report
def testFASTAFile():
  """
  @tests: FASTAFile
   FASTAFile.__init__
   FASTAFile.load
   FASTAFile.store
  """
  ff = pyopenms.FASTAFile()

  assert pyopenms.FASTAFile().load is not None
  assert pyopenms.FASTAFile().store is not None


@report
def testFASTAEntry():
  """
  @tests: FASTAEntry
   FASTAEntry.__init__
  """
  ff = pyopenms.FASTAEntry()


@report
def testTransformationXMLFile():
  """
  @tests: TransformationXMLFile
   TransformationXMLFile.__init__
   TransformationXMLFile.load
   TransformationXMLFile.store
  """
  fh = pyopenms.TransformationXMLFile()
  td = pyopenms.TransformationDescription()
  fh.store("test.transformationXML", td)
  fh.load("test.transformationXML", td, True)
  assert td.getDataPoints() == []


@report
def testParamXMLFile():
  """
  @tests: ParamXMLFile
   ParamXMLFile.__init__
   ParamXMLFile.load
   ParamXMLFile.store
  """
  fh = pyopenms.ParamXMLFile()
  p = pyopenms.Param()
  fh.store("test.ini", p)
  fh.load("test.ini", p)


@report
def testTransitionTSVFile():
  """
  @tests: TransitionTSVFile
   TransitionTSVFile.__init__
   TransitionTSVFile.convertTargetedExperimentToTSV
   TransitionTSVFile.convertTSVToTargetedExperiment
   TransitionTSVFile.validateTargetedExperiment
  """
  ff = pyopenms.TransitionTSVFile()

  assert pyopenms.TransitionTSVFile().convertTargetedExperimentToTSV is not None
  assert pyopenms.TransitionTSVFile().convertTSVToTargetedExperiment is not None
  assert pyopenms.TransitionTSVFile().validateTargetedExperiment is not None


@report
def testSwathFile():
  """
  @tests: SwathFile
   SwathFile.__init__
  """
  fh = pyopenms.SwathFile()


@report
def testIBSpectraFile():
  """
  @tests: IBSpectraFile
   IBSpectraFile.__init__
   IBSpectraFile.store
  """
  fh = pyopenms.IBSpectraFile()
  cmap = pyopenms.ConsensusMap()
  correctError = False
  try:
    fh.store(pyopenms.String("test.ibspectra.file"), cmap)
    raise AssertionError("Should have raised RuntimeError")
  except RuntimeError:
    correctError = True

  assert correctError


@report
def testCVMappings():
  """
  @tests: CVMappings
   CVMappings.__init__
  """
  val = pyopenms.CVMappings()


@report
def testCVMappingFile():
  """
  @tests: CVMappingFile
   CVMappingFile.__init__
   CVMappingFile.load
  """
  val = pyopenms.CVMappingFile()

  assert pyopenms.CVMappingFile().load


@report
def testControlledVocabulary():
  """
  @tests: ControlledVocabulary
   ControlledVocabulary.__init__
   ControlledVocabulary.loadFromOBO
  """
  val = pyopenms.ControlledVocabulary()

  assert pyopenms.ControlledVocabulary().loadFromOBO


if __name__ == "__main__":
  # Run all tests in this module
  import sys
  module = sys.modules[__name__]
    
  test_functions = [getattr(module, name) for name in dir(module) 
           if name.startswith('test') and callable(getattr(module, name))]
    
  print(f"Running {len(test_functions)} file I/O tests...")
  for test_func in test_functions:
    try:
      test_func()
      print(f"✅ {test_func.__name__}")
    except Exception as e:
      print(f"❌ {test_func.__name__}: {e}")
