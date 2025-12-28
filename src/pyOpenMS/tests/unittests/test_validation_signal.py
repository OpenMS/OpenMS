"""
Test file for validation and signal processing functionality.

Part of test000.py split (Issue #8567).
Contains tests for:
- Signal-to-noise estimation
- Validation (semantic validator, checksum types)
- ConsensusID algorithms
- Internal calibration
- File handling (MzML/MzXML/MzData, PeakFileOptions)
- Percolator output

Tests: 14 total
"""

import pyopenms as oms
from test_helpers import report, _testParam, _testProgressLogger


@report
def testChecksumType():
  """
  @tests: ChecksumType
   ChecksumType.MD5
   ChecksumType.SHA1
   ChecksumType.SIZE_OF_CHECKSUMTYPE
   ChecksumType.UNKNOWN_CHECKSUM
  """
  assert isinstance(oms.ChecksumType.MD5, int)
  assert isinstance(oms.ChecksumType.SHA1, int)
  assert isinstance(oms.ChecksumType.SIZE_OF_CHECKSUMTYPE, int)
  assert isinstance(oms.ChecksumType.UNKNOWN_CHECKSUM, int)


@report
def testSignalToNoiseEstimatorMedian():
  """
  @tests: SignalToNoiseEstimatorMedian
   SignalToNoiseEstimatorMedian.__init__
  """
  f = oms.SignalToNoiseEstimatorMedian()
  assert f.init is not None
  assert f.getSignalToNoise is not None


@report
def testSignalToNoiseEstimatorMedianChrom():
  """
  @tests: SignalToNoiseEstimatorMedianChrom
   SignalToNoiseEstimatorMedianChrom.__init__
  """
  f = oms.SignalToNoiseEstimatorMedianChrom()
  assert f.init is not None
  assert f.getSignalToNoise is not None


@report
def testSemanticValidator():
  """
  @tests: SemanticValidator
   SemanticValidator.__init__
  """
  m = oms.CVMappings()
  cv = oms.ControlledVocabulary()

  val = oms.SemanticValidator(m, cv)

  assert val.validate is not None
  assert val.setCheckTermValueTypes is not None
  assert val.setCheckUnits is not None


@report
def testInternalCalibration():
  """
  @tests: InternalCalibration
   InternalCalibration.__init__
   InternalCalibration.calibrateMapGlobally
   InternalCalibration.calibrateMapSpectrumwise
  """
  ff = oms.InternalCalibration()

  assert oms.InternalCalibration().fillCalibrants is not None
  assert oms.InternalCalibration().getCalibrationPoints is not None
  assert oms.InternalCalibration().calibrate is not None


@report
def testMxxxFile():
  """
  @tests: MzDataFile
   MzDataFile.__init__
   MzDataFile.endProgress
   MzDataFile.getLogType
   MzDataFile.load
   MzDataFile.setLogType
   MzDataFile.setProgress
   MzDataFile.startProgress
   MzDataFile.store
   MzDataFile.getOptions
   MzDataFile.setOptions

   MzMLFile.__init__
   MzMLFile.endProgress
   MzMLFile.getLogType
   MzMLFile.load
   MzMLFile.setLogType
   MzMLFile.setProgress
   MzMLFile.startProgress
   MzMLFile.store
   MzMLFile.getOptions
   MzMLFile.setOptions

   MzXMLFile.getOptions
   MzXMLFile.setOptions
   MzXMLFile.__init__
   MzXMLFile.endProgress
   MzXMLFile.getLogType
   MzXMLFile.load
   MzXMLFile.setLogType
   MzXMLFile.setProgress
   MzXMLFile.startProgress
   MzXMLFile.store

  """
  mse = oms.MSExperiment()
  s = oms.MSSpectrum()
  mse.addSpectrum(s)

  fh = oms.MzDataFile()
  _testProgressLogger(fh)
  fh.store("test.mzData", mse)
  fh.load("test.mzData", mse)

  fh.setOptions(fh.getOptions())

  fh = oms.MzMLFile()
  _testProgressLogger(fh)
  fh.store("test.mzML", mse)
  fh.load("test.mzML", mse)
  fh.setOptions(fh.getOptions())

  myStr = oms.String()
  fh.storeBuffer(myStr, mse)
  assert len(myStr.toString()) == 5269
  mse2 = oms.MSExperiment()
  fh.loadBuffer(bytes(myStr), mse2)
  assert mse2 == mse
  assert mse2.size() == 1

  fh = oms.MzXMLFile()
  _testProgressLogger(fh)
  fh.store("test.mzXML", mse)
  fh.load("test.mzXML", mse)
  fh.setOptions(fh.getOptions())


@report
def testPeakFileOptions():
  """
  @tests: PeakFileOptions
   PeakFileOptions.__init__
   PeakFileOptions.addMSLevel
   PeakFileOptions.clearMSLevels
   PeakFileOptions.containsMSLevel
   PeakFileOptions.getCompression
   PeakFileOptions.getMSLevels
   PeakFileOptions.getMetadataOnly
   PeakFileOptions.getWriteSupplementalData
   PeakFileOptions.hasMSLevels
   PeakFileOptions.setCompression
   PeakFileOptions.setMSLevels
   PeakFileOptions.setMetadataOnly
   PeakFileOptions.setWriteSupplementalData
  """

  pfo = oms.PeakFileOptions()
  pfo.addMSLevel
  pfo.clearMSLevels()
  pfo.containsMSLevel(1)
  pfo.getCompression()
  pfo.getMSLevels()
  pfo.getMetadataOnly()
  pfo.getWriteSupplementalData()
  pfo.hasMSLevels()
  pfo.setCompression
  pfo.setMSLevels
  pfo.setMetadataOnly
  pfo.setWriteSupplementalData


@report
def testConsensusIDAlgorithmAverage():
  algo = oms.ConsensusIDAlgorithmAverage()
  assert algo.apply


@report
def testConsensusIDAlgorithmBest():
  algo = oms.ConsensusIDAlgorithmBest()
  assert algo.apply


@report
def testConsensusIDAlgorithmPEPIons():
  algo = oms.ConsensusIDAlgorithmPEPIons()
  assert algo.apply


@report
def testConsensusIDAlgorithmPEPMatrix():
  algo = oms.ConsensusIDAlgorithmPEPMatrix()
  assert algo.apply


@report
def testConsensusIDAlgorithmRanks():
  algo = oms.ConsensusIDAlgorithmRanks()
  assert algo.apply


@report
def testConsensusIDAlgorithmWorst():
  algo = oms.ConsensusIDAlgorithmWorst()
  assert algo.apply


@report
def testPercolatorOutfile():
  e = oms.PercolatorOutfile()
  assert e


if __name__ == '__main__':
  print("Running 14 validation and signal processing tests...")
    
  tests = [
    testChecksumType,
    testSignalToNoiseEstimatorMedian,
    testSignalToNoiseEstimatorMedianChrom,
    testSemanticValidator,
    testInternalCalibration,
    testMxxxFile,
    testPeakFileOptions,
    testConsensusIDAlgorithmAverage,
    testConsensusIDAlgorithmBest,
    testConsensusIDAlgorithmPEPIons,
    testConsensusIDAlgorithmPEPMatrix,
    testConsensusIDAlgorithmRanks,
    testConsensusIDAlgorithmWorst,
    testPercolatorOutfile,
  ]
    
  passed = 0
  failed = 0
    
  for test in tests:
    try:
      test()
      passed += 1
      print(f"✅ {test.__name__}")
    except Exception as e:
      failed += 1
      print(f"❌ {test.__name__}: {e}")
    
  print(f"\nResults: {passed} passed, {failed} failed")

