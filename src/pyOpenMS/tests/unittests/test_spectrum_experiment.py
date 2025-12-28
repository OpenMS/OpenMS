#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Spectrum and Experiment tests extracted from test000.py
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
  _testParam,
  _testStrOutput,
  extraction_performance_test
)

@report
def testSpectrumAlignment():
  """
  @tests: SpectrumAlignment
    SpectrumAlignment.__init__
    SpectrumAlignment.getSpectrumAlignment
  """
  # test existence of some methods
  pyopenms.SpectrumAlignment
  pyopenms.SpectrumAlignment.__init__
  pyopenms.SpectrumAlignment.getDefaults
  pyopenms.SpectrumAlignment.getParameters
  pyopenms.SpectrumAlignment.setParameters

  spec = pyopenms.MSSpectrum()
  p = pyopenms.Peak1D()
  p.setMZ(1000.0)
  p.setIntensity(200.0)
  spec.push_back(p)
  p.setMZ(2000.0)
  p.setIntensity(200.0)
  spec.push_back(p)

  rich_spec = pyopenms.MSSpectrum()
  p = pyopenms.Peak1D()
  p.setMZ(1000.001)
  p.setIntensity(200.0)
  rich_spec.push_back(p)
  p.setMZ(2000.001)
  p.setIntensity(200.0)
  rich_spec.push_back(p)
  p.setMZ(3000.001)
  p.setIntensity(200.0)
  rich_spec.push_back(p)

  aligner = pyopenms.SpectrumAlignment()
  result = []

  aligner.getSpectrumAlignment(result, spec, spec)
  assert result == [ (0,0), (1,1) ], result
  aligner.getSpectrumAlignment(result, rich_spec, spec)
  assert result == [ (0,0), (1,1) ], result
  aligner.getSpectrumAlignment(result, spec, rich_spec)
  assert result == [ (0,0), (1,1) ], result
  aligner.getSpectrumAlignment(result, rich_spec, rich_spec)
  assert result == [ (0,0), (1,1), (2,2) ], result

  aligner = pyopenms.SpectrumAlignmentScore()
  assert isinstance(aligner(spec), float)
  assert isinstance(aligner(rich_spec), float)

  assert isinstance(aligner(spec, rich_spec), float)
  assert isinstance(aligner(rich_spec, spec), float)

  assert isinstance(aligner(spec, spec), float)
  assert isinstance(aligner(rich_spec, rich_spec), float)

@report
def testChromatogramTools():
  """
  @tests: ChromatogramTools
   ChromatogramTools.__init__
   ChromatogramTools.convertChromatogramsToSpectra
   ChromatogramTools.convertSpectraToChromatograms
  """
  pyopenms.ChromatogramTools()
  assert pyopenms.ChromatogramTools.convertChromatogramsToSpectra is not None
  assert pyopenms.ChromatogramTools.convertSpectraToChromatograms is not None


@report
def testExperimentalSettings():
  """
  @tests: ExperimentalSettings
   ExperimentalSettings.__init__
  """
  pyopenms.ExperimentalSettings()

@report
def testMSExperiment():
  """
  @tests: MSExperiment
   MSExperiment.__init__
   MSExperiment.getLoadedFilePath
   MSExperiment.getMaxMZ
   MSExperiment.getMaxRT
   MSExperiment.getMetaValue
   MSExperiment.getMinMZ
   MSExperiment.getMinRT
   MSExperiment.push_back
   MSExperiment.setLoadedFilePath
   MSExperiment.setMetaValue
   MSExperiment.size
   MSExperiment.sortSpectra
   MSExperiment.updateRanges
   MSExperiment.__eq__
   MSExperiment.__ge__
   MSExperiment.__getitem__
   MSExperiment.__gt__
   MSExperiment.__iter__
   MSExperiment.__le__
   MSExperiment.__lt__
   MSExperiment.__ne__
   MSExperiment.clearMetaInfo
   MSExperiment.getKeys
   MSExperiment.isMetaEmpty
   MSExperiment.metaValueExists
   MSExperiment.removeMetaValue
   MSExperiment.getSize
   MSExperiment.isSorted
   MSExperiment.get2DPeakDataLong
   MSExperiment.get_df
   MSExperiment.get_massql_df
  """
  mse = pyopenms.MSExperiment()
  mse_ = copy.copy(mse)
  assert mse_ == mse
  mse_ = copy.deepcopy(mse)
  assert mse_ == mse
  mse_ = pyopenms.MSExperiment(mse)
  assert mse_ == mse

  _testMetaInfoInterface(mse)

  mse.setLoadedFilePath("")
  assert mse.size() == 0

  mse.getIdentifier()
  mse.getLoadedFileType()
  mse.getLoadedFilePath()

  # We need to add a feature to the map before updateRanges otherwise the getMin and getMax throw an error.
  spec = pyopenms.MSSpectrum()
  data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
  data_i = np.array( [50.0, 80.0] ).astype(np.float32)
  spec.set_peaks( [data_mz,data_i] )

  mse.addSpectrum(spec)
  assert mse.size() == 1

  mse.updateRanges()
  mse.sortSpectra(True)

  assert isinstance(mse.getMaxRT(), float)
  assert isinstance(mse.getMinRT(), float)
  assert isinstance(mse.getMaxMZ(), float)
  assert isinstance(mse.getMinMZ(), float)
  _testStrOutput(mse.getLoadedFilePath())
  assert isinstance(mse.getMinIntensity(), float)
  assert isinstance(mse.getMaxIntensity(), float)

  assert mse[0] is not None

  mse.updateRanges()
  rt, mz, inty = mse.get2DPeakDataLong(mse.getMinRT(), mse.getMaxRT(), mse.getMinMZ(), mse.getMaxMZ(), 1)
  assert rt.shape[0] == 2
  assert mz.shape[0] == 2
  assert inty.shape[0] == 2

  assert isinstance(list(mse), list)

  assert mse == mse
  assert not mse != mse

  assert mse.getSize() >= 0
  assert int(mse.isSorted()) in (0,1)

  mse2 = copy.copy(mse)

  assert mse.getSize() == mse2.getSize()
  assert mse2 == mse

  exp = pyopenms.MSExperiment()

  for i in range(5):
    s = pyopenms.MSSpectrum()
    s.setRT(i)
    s.setMSLevel(1 if i % 2 == 0 else 2)

    for mz in (500, 600):
      p = pyopenms.Peak1D()
      p.setMZ(mz + i)
      p.setIntensity(i + 10)
      s.push_back(p)

    exp.addSpectrum(s)

  assert exp.get_df().shape == (5, 4)
  assert exp.get_df(ms_levels=[1]).shape == (3, 4)
  assert exp.get_df(ms_levels=[2]).shape == (2, 4)

  assert exp.get_df(long_format=True).shape == (10, 4)
  assert exp.get_df(long_format=True, ms_levels=[1]).shape == (6, 4)
  assert exp.get_df(long_format=True, ms_levels=[2]).shape == (4, 4)

  pyopenms.MzMLFile().load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BSA1_F1.mzML'), exp)

  ms1_df, ms2_df = exp.get_massql_df()
  assert ms1_df.shape == (140055, 7)
  assert np.allclose(ms2_df.head(), pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BSA1_F1_MS2_MassQL.tsv'), sep='\t'))

  pyopenms.MzMLFile().load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BSA1_F1_ION.mzML'), exp)
  df = exp.get_ion_df()
  assert np.allclose(df.head(), pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BSA1_F1_MS1_ION.tsv'), sep='\t'))

  ms1_df, ms2_df = exp.get_massql_df(ion_mobility=True)
  assert ms1_df.shape == (332620, 8)
  assert np.allclose(ms1_df.head(), pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BSA1_F1_MS1_MassQL_ION.tsv'), sep='\t'))

  #####################################################################################
  # test fast aggregation and XIC extraction using ranges
  pyopenms.MzMLFile().load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BSA1_F1.mzML'), exp)    
  exp.updateRanges()

  ############################################################################
  # Uncomment to run performance tests
  # print("\nStarting performance tests MS1...")
  # aggregate_time, xic_time = extraction_performance_test(exp, 1) # Run performance tests on MS level 1   
  # print("\nPerformance Results:")
  # print(f"aggregateFromMatrix time for 1M ranges: {aggregate_time:.2f} seconds")
  # print(f"extractXICsFromMatrix time for 1M ranges: {xic_time:.2f} seconds")
  # print(f"Ratio (XIC/aggregate): {xic_time/aggregate_time:.2f}")

  # print("\nStarting performance tests MS2...")
  # aggregate_time, xic_time = extraction_performance_test(exp, 2) # Run performance tests on MS level 1   
  # print("\nPerformance Results:")
  # print(f"aggregateFromMatrix time for 1M ranges: {aggregate_time:.2f} seconds")
  # print(f"extractXICsFromMatrix time for 1M ranges: {xic_time:.2f} seconds")
  # print(f"Ratio (XIC/aggregate): {xic_time/aggregate_time:.2f}")
  # assert false # needed to output the results

  # eluting peptide feature at these coordinates
  rt_min = 1730.0
  rt_max = 1770.0
  iso1_mz = 443.711 # monoisotopic peak
  iso2_mz = 444.212 # first isotopic peak
  iso3_mz = 444.713 # second isotopic peak

  no_iso1_mz = 444.000 # no peaks in this area
  no_iso2_mz = 444.403 # some noise peaks in this area


  # Create ranges matrix with structure:
  # [[mz_min, mz_max, rt_min, rt_max], ...]
  ranges_matrix = pyopenms.MatrixDouble.fromNdArray(
     np.array([
    # Expected isotope peaks
    [iso1_mz - 0.1, iso1_mz + 0.1, rt_min, rt_max],
    [iso2_mz - 0.1, iso2_mz + 0.1, rt_min, rt_max],
    [iso3_mz - 0.1, iso3_mz + 0.1, rt_min, rt_max],
    # Regions where we don't expect peaks
    [no_iso1_mz - 0.1, no_iso1_mz + 0.1, rt_min, rt_max],
    [no_iso2_mz - 0.1, no_iso2_mz + 0.1, rt_min, rt_max]
  ])
  )
    
  # Print ranges matrix for verification
  print("Ranges Matrix:")
  print (ranges_matrix.get_matrix())

  agg_result = exp.aggregateFromMatrix(ranges_matrix, 1, b"sum")    
  agg_result_array = np.array(agg_result)    # Convert result to numpy array for easier testing

  print("\nAggregation Results:")
  print(agg_result)

  # Basic shape checks
  assert len(agg_result_array) == ranges_matrix.rows(), f"Expected {ranges_matrix.rows()} results, got {len(agg_result_array)}"
    
  # Check that regions with expected peaks have higher values
  isotope_intensities = agg_result_array[:3]  # First three rows are isotope peaks
  no_peak_intensities = agg_result_array[3:]  # Last two rows are regions without peaks
    
  print(isotope_intensities)
  print(no_peak_intensities)
 
  print("\nIsotope Peak Arrays and Their Sums:")
  isotope_sums = []
  for i, intensity_array in enumerate(isotope_intensities):
    array_sum = np.sum(intensity_array)
    isotope_sums.append(array_sum)
    print(f"Isotope {i+1} array: {intensity_array}")
    print(f"Isotope {i+1} sum: {array_sum}")
    
  print("\nNo-Peak Region Arrays and Their Sums:")
  no_peak_sums = []
  for i, intensity_array in enumerate(no_peak_intensities):
    array_sum = np.sum(intensity_array)
    no_peak_sums.append(array_sum)
    print(f"No-peak region {i+1} array: {intensity_array}")
    print(f"No-peak region {i+1} sum: {array_sum}")

  EXPECTED_ISO_SUMS = [24680058.50, 11043987.94, 3141677.76] 
  EXPECTED_NO_PEAK_SUMS = [0.0, 12322.06]  
    
  # Test each isotope array sum
  for i, (actual_sum, expected_sum) in enumerate(zip(isotope_sums, EXPECTED_ISO_SUMS, strict=True)):
    assert np.isclose(actual_sum, expected_sum, rtol=1e-5), \
      f"Expected isotope {i+1} sum {expected_sum}, got {actual_sum}"
    
  # Test each no-peak array sum
  for i, (actual_sum, expected_sum) in enumerate(zip(no_peak_sums, EXPECTED_NO_PEAK_SUMS, strict=True)):
    assert np.isclose(actual_sum, expected_sum, rtol=1e-5), \
      f"Expected no-peak region {i+1} sum {expected_sum}, got {actual_sum}"


  xic_result = exp.extractXICsFromMatrix(ranges_matrix, 1, b"sum")
  print("\nXIC Results:")
  for chrom in xic_result:
    print(f"m/z: {chrom.getProduct().getMZ()}")
    print(f"Size: {chrom.size()}")

  xic_details = []
  for i, chrom in enumerate(xic_result):
    # Get basic XIC information
    mz = chrom.getProduct().getMZ()
    size = chrom.size()
        
    # Get intensity values from the chromatogram
    intensities = [point.getIntensity() for point in chrom]
    rt_values = [point.getRT() for point in chrom]
    total_intensity = sum(intensities)
        
    details = {
      'mz': mz,
      'size': size,
      'total_intensity': total_intensity,
      'rt_range': (min(rt_values), max(rt_values)) if rt_values else (None, None)
    }
    xic_details.append(details)
        
    print(f"\nXIC {i+1}:")
    print(f"m/z: {mz}")
    print(f"Number of points: {size}")
    print(f"Total intensity: {total_intensity}")
    print(f"RT range: {details['rt_range']}")

  # XIC expectations
  EXPECTED_XIC_SIZES = [24, 24, 24, 24, 24]  # Expected number of points in each XIC
  EXPECTED_XIC_TOTAL_INTENSITIES = [24680058.50, 11043987.94, 3141677.76, 0.0, 12322.06]  # Expected total intensity for each XIC
    
  # Test XIC results
  for i, (details, exp_size, exp_total) in enumerate(zip(
      xic_details, EXPECTED_XIC_SIZES, EXPECTED_XIC_TOTAL_INTENSITIES, strict=True)):
        
    assert details['size'] == exp_size, \
      f"XIC {i+1}: Expected {exp_size} points, got {details['size']}"
        
    assert np.isclose(details['total_intensity'], exp_total, rtol=1e-5), \
      f"XIC {i+1}: Expected total intensity {exp_total}, got {details['total_intensity']}"
                
    # Verify RT range falls within expected bounds
    if details['rt_range'][0] is not None:
      assert rt_min <= details['rt_range'][0] <= rt_max, \
        f"XIC {i+1}: Start RT {details['rt_range'][0]} outside expected range [{rt_min}, {rt_max}]"
      assert rt_min <= details['rt_range'][1] <= rt_max, \
        f"XIC {i+1}: End RT {details['rt_range'][1]} outside expected range [{rt_min}, {rt_max}]"

  # Test __repr__ and __str__ methods for MSExperiment
  exp_repr = pyopenms.MSExperiment()
  s1 = pyopenms.MSSpectrum()
  s1.setRT(100.0)
  s1.setMSLevel(1)
  s1.set_peaks(([100.0, 200.0], [1000.0, 2000.0]))
  s2 = pyopenms.MSSpectrum()
  s2.setRT(200.0)
  s2.setMSLevel(2)
  s2.set_peaks(([150.0, 250.0], [500.0, 1500.0]))
  exp_repr.addSpectrum(s1)
  exp_repr.addSpectrum(s2)
  exp_repr.updateRanges()

  repr_str = repr(exp_repr)
  assert "MSExperiment(" in repr_str
  assert "num_spectra=" in repr_str
  assert "num_chromatograms=" in repr_str
  str_str = str(exp_repr)
  assert str_str == repr_str


def testSpectrumSetting(s=None):
  """
  @tests: SpectrumSettings
   SpectrumSettings.SpectrumType
   SpectrumSettings.__init__
   SpectrumSettings.getAcquisitionInfo
   SpectrumSettings.getComment
   SpectrumSettings.getDataProcessing
   SpectrumSettings.getInstrumentSettings
   SpectrumSettings.getNativeID
   SpectrumSettings.getPrecursors
   SpectrumSettings.getProducts
   SpectrumSettings.getSourceFile
   SpectrumSettings.getType
   SpectrumSettings.setAcquisitionInfo
   SpectrumSettings.setComment
   SpectrumSettings.setDataProcessing
   SpectrumSettings.setInstrumentSettings
   SpectrumSettings.setNativeID
   SpectrumSettings.setPrecursors
   SpectrumSettings.setProducts
   SpectrumSettings.setSourceFile
   SpectrumSettings.setType
   SpectrumSettings.unify
  """
  if s is None:
    s = pyopenms.SpectrumSettings()

  assert s.getType() in [ pyopenms.SpectrumSettings.SpectrumType.UNKNOWN,
                 pyopenms.SpectrumSettings.SpectrumType.CENTROID,
                 pyopenms.SpectrumSettings.SpectrumType.PROFILE]

  assert isinstance(s.getAcquisitionInfo(), pyopenms.AcquisitionInfo)
  assert isinstance(s.getInstrumentSettings(), pyopenms.InstrumentSettings)
  assert isinstance(s.getSourceFile(), pyopenms.SourceFile)
  assert isinstance(s.getDataProcessing(), list)

  s.setAcquisitionInfo(s.getAcquisitionInfo())
  s.setInstrumentSettings(s.getInstrumentSettings())
  s.setSourceFile(s.getSourceFile())
  s.setDataProcessing(s.getDataProcessing())
  s.setComment(s.getComment())
  s.setPrecursors(s.getPrecursors())
  s.setProducts(s.getProducts())
  s.setType(s.getType())
  s.setNativeID(s.getNativeID())
  s.setType(s.getType())
  if isinstance(s, pyopenms.SpectrumSettings):
    s.unify(s)

  # Test getAllNamesOf method
  spectrum_type_names = pyopenms.SpectrumSettings.getAllNamesOfSpectrumType()
  assert len(spectrum_type_names) == pyopenms.SpectrumSettings.SpectrumType.SIZE_OF_SPECTRUMTYPE
  assert spectrum_type_names[pyopenms.SpectrumSettings.SpectrumType.CENTROID].decode() == "Centroid"
  assert spectrum_type_names[pyopenms.SpectrumSettings.SpectrumType.PROFILE].decode() == "Profile"


@report
def testMSSpectrum():
  """
  @tests: MSSpectrum
   MSSpectrum.__init__
   MSSpectrum.clear
   MSSpectrum.clearMetaInfo
   MSSpectrum.findNearest
   MSSpectrum.getAcquisitionInfo
   MSSpectrum.getComment
   MSSpectrum.getDataProcessing
   MSSpectrum.getInstrumentSettings
   MSSpectrum.getKeys
   MSSpectrum.getMSLevel
   MSSpectrum.getMetaValue
   MSSpectrum.getName
   MSSpectrum.getNativeID
   MSSpectrum.getPrecursors
   MSSpectrum.getProducts
   MSSpectrum.getRT
   MSSpectrum.getSourceFile
   MSSpectrum.getType
   MSSpectrum.get_peaks
   MSSpectrum.intensityInRange
   MSSpectrum.isMetaEmpty
   MSSpectrum.isSorted
   MSSpectrum.metaValueExists
   MSSpectrum.push_back
   MSSpectrum.removeMetaValue
   MSSpectrum.setAcquisitionInfo
   MSSpectrum.setComment
   MSSpectrum.setDataProcessing
   MSSpectrum.setInstrumentSettings
   MSSpectrum.setMSLevel
   MSSpectrum.setMetaValue
   MSSpectrum.setName
   MSSpectrum.setNativeID
   MSSpectrum.setPeptideIdentifications
   MSSpectrum.setPrecursors
   MSSpectrum.setProducts
   MSSpectrum.setRT
   MSSpectrum.setSourceFile
   MSSpectrum.setType
   MSSpectrum.set_peaks
   MSSpectrum.size
   MSSpectrum.unify
   MSSpectrum.updateRanges
   MSSpectrum.__eq__
   MSSpectrum.__ge__
   MSSpectrum.__getitem__
   MSSpectrum.__gt__
   MSSpectrum.__le__
   MSSpectrum.__lt__
   MSSpectrum.__ne__
   """
  spec = pyopenms.MSSpectrum()
  spec_ = copy.copy(spec)
  assert spec_ == spec
  spec_ = copy.deepcopy(spec)
  assert spec_ == spec
  spec_ = pyopenms.MSSpectrum(spec)
  assert spec_ == spec

  _testMetaInfoInterface(spec)

  testSpectrumSetting(spec)

  spec.setRT(3.0)
  assert spec.getRT() == 3.0
  spec.setMSLevel(2)
  assert spec.getMSLevel() == 2
  spec.setName("spec")
  assert spec.getName() == "spec"

  p = pyopenms.Peak1D()
  p.setMZ(1000.0)
  p.setIntensity(200.0)
  spec.push_back(p)

  assert spec.size() == 1
  assert spec[0] == p

  spec.updateRanges()
  assert isinstance(spec.findNearest(0.0), int)

  assert isinstance(spec.getMinMZ(), float)
  assert isinstance(spec.getMaxMZ(), float)
  assert isinstance(spec.getMinIntensity(), float)
  assert isinstance(spec.getMaxIntensity(), float)

  assert spec == spec
  assert not spec != spec

  mz, ii = spec.get_peaks()
  assert len(mz) == len(ii)
  assert len(mz) == 1

  spec.set_peaks((mz, ii))
  mz0, ii0 = spec.get_peaks()
  assert mz0 == mz
  assert ii0 == ii

  assert int(spec.isSorted()) in (0,1)

  spec.clear(False)
  p = pyopenms.Peak1D()
  p.setMZ(1000.0)
  p.setIntensity(200.0)
  spec.push_back(p)
  p = pyopenms.Peak1D()
  p.setMZ(2000.0)
  p.setIntensity(400.0)
  spec.push_back(p)

  mz, ii = spec.get_peaks()
  assert spec[0].getMZ() == 1000.0
  assert spec[1].getMZ() == 2000.0
  assert spec[0].getIntensity() == 200.0
  assert spec[1].getIntensity() == 400.0
  assert mz[0] == 1000.0
  assert mz[1] == 2000.0
  assert ii[0] == 200.0
  assert ii[1] == 400.0

  spec.setMSLevel(2)
  spec.setNativeID('scan=1')
  prec = pyopenms.Precursor()
  prec.setMZ(100.0)
  prec.setCharge(1)
  spec.setPrecursors([prec])
  spec.setMetaValue('total ion current', 600)

  data = np.array( [5, 8] ).astype(np.float32)
  f_da = [ pyopenms.FloatDataArray() ]
  f_da[0].set_data(data)
  f_da[0].setName("Ion Mobility")
  spec.setFloatDataArrays( f_da )
  spec.setDriftTimeUnit( pyopenms.DriftTimeUnit.MILLISECOND )

  s_da = pyopenms.StringDataArray()
  for s in ['b3+', 'y4+']:
    s_da.push_back(s)
  s_da.setName("IonNames")
  spec.setStringDataArrays([s_da])

  df = spec.get_df()
  assert df.shape == (2, 11)
  assert df.loc[0, 'mz'] == 1000.0
  assert df.loc[0, 'intensity'] == 200.0
  assert df.loc[0, 'ion_mobility'] == 5.0
  assert df.loc[0, 'ion_mobility_unit'] == 'ms'
  assert df.loc[0, 'ms_level'] == 2
  assert df.loc[0, 'precursor_mz'] == 100.0
  assert df.loc[0, 'precursor_charge'] == 1
  assert df.loc[0, 'native_id'] == 'scan=1'
  assert df.loc[0, 'ion_annotation'] == 'b3+'
  assert df.loc[0, 'total ion current'] == 600

  spec.clear(False)
  data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
  data_i = np.array( [50.0, 80.0] ).astype(np.float32)
  spec.set_peaks( [data_mz,data_i] )

  mz, ii = spec.get_peaks()
  assert spec[0].getMZ() == 5.0
  assert spec[1].getMZ() == 8.0
  assert spec[0].getIntensity() == 50.0
  assert spec[1].getIntensity() == 80.0
  assert mz[0] == 5.0
  assert mz[1] == 8.0
  assert ii[0] == 50.0
  assert ii[1] == 80.0

  # Fast
  spec.clear(False)
  data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
  data_i = np.array( [50.0, 80.0] ).astype(np.float64)
  spec.set_peaks( [data_mz,data_i] )

  mz, ii = spec.get_peaks()
  assert spec[0].getMZ() == 5.0
  assert spec[1].getMZ() == 8.0
  assert spec[0].getIntensity() == 50.0
  assert spec[1].getIntensity() == 80.0
  assert mz[0] == 5.0
  assert mz[1] == 8.0
  assert ii[0] == 50.0
  assert ii[1] == 80.0

  # Slow
  spec.clear(False)
  data_mz = np.array( [5.0, 8.0] ).astype(np.float32)
  data_i = np.array( [50.0, 80.0] ).astype(np.float32)
  spec.set_peaks( [data_mz,data_i] )

  mz, ii = spec.get_peaks()
  assert spec[0].getMZ() == 5.0
  assert spec[1].getMZ() == 8.0
  assert spec[0].getIntensity() == 50.0
  assert spec[1].getIntensity() == 80.0
  assert mz[0] == 5.0
  assert mz[1] == 8.0
  assert ii[0] == 50.0
  assert ii[1] == 80.0

  ###################################
  # get data arrays
  ###################################
  assert len(spec.getStringDataArrays()) == 0
  string_da = [ pyopenms.StringDataArray() ]
  string_da[0].push_back("hello")
  string_da[0].push_back("world")
  string_da[0].setName("greetings")
  string_da.append( pyopenms.StringDataArray() )
  string_da[1].push_back("other")
  spec.setStringDataArrays( string_da )
  assert len(spec.getStringDataArrays()) == 2
  assert spec.getStringDataArrays()[0][0] == b"hello"
  assert spec.getStringDataArrays()[1][0] == b"other"
  assert spec.getStringDataArrays()[0] == spec.getStringDataArrays()[0] # test __eq__
  assert spec.getStringDataArrays()[0] != spec.getStringDataArrays()[1] # test __ne__


  spec = pyopenms.MSSpectrum()
  assert len(spec.getIntegerDataArrays()) == 0
  int_da = [ pyopenms.IntegerDataArray() ]
  int_da[0].push_back(5)
  int_da[0].push_back(6)
  int_da[0].setName("test")
  int_da.append( pyopenms.IntegerDataArray() )
  int_da[1].push_back(8)
  spec.setIntegerDataArrays( int_da )
  assert len(spec.getIntegerDataArrays()) == 2
  assert spec.getIntegerDataArrays()[0][0] == 5
  assert spec.getIntegerDataArrays()[1][0] == 8
  assert spec.getIntegerDataArrays()[0] == spec.getIntegerDataArrays()[0] # test __eq__
  assert spec.getIntegerDataArrays()[0] != spec.getIntegerDataArrays()[1] # test __ne__

  spec = pyopenms.MSSpectrum()
  data = np.array( [5, 8, 42] ).astype(np.intc)
  int_da = [ pyopenms.IntegerDataArray() ]
  int_da[0].set_data(data)
  spec.setIntegerDataArrays( int_da )
  assert len(spec.getIntegerDataArrays()) == 1
  assert spec.getIntegerDataArrays()[0][0] == 5
  assert spec.getIntegerDataArrays()[0][2] == 42
  assert len(int_da[0].get_data() ) == 3

  spec = pyopenms.MSSpectrum()
  assert len(spec.getFloatDataArrays()) == 0
  f_da = [ pyopenms.FloatDataArray() ]
  f_da[0].push_back(5.0)
  f_da[0].push_back(6.0)
  f_da[0].setName("test")
  f_da.append( pyopenms.FloatDataArray() )
  f_da[1].push_back(8.0)
  spec.setFloatDataArrays( f_da )
  assert len(spec.getFloatDataArrays()) == 2
  assert spec.getFloatDataArrays()[0][0] == 5.0
  assert spec.getFloatDataArrays()[1][0] == 8.0
  assert spec.getFloatDataArrays()[0] == spec.getFloatDataArrays()[0] # test __eq__
  assert spec.getFloatDataArrays()[0] != spec.getFloatDataArrays()[1] # test __ne__

  spec = pyopenms.MSSpectrum()
  data = np.array( [5, 8, 42] ).astype(np.float32)
  f_da = [ pyopenms.FloatDataArray() ]
  f_da[0].set_data(data)
  spec.setFloatDataArrays( f_da )
  assert len(spec.getFloatDataArrays()) == 1
  assert spec.getFloatDataArrays()[0][0] == 5.0
  assert spec.getFloatDataArrays()[0][2] == 42.0
  assert len(f_da[0].get_data() ) == 3

  spec = pyopenms.MSSpectrum()
  dfunit = spec.getDriftTimeUnit()
  assert pyopenms.DriftTimeUnit().getMapping()[dfunit] == "NONE"
  assert dfunit == pyopenms.DriftTimeUnit.NONE
  assert spec.getDriftTimeUnitAsString() == '<NONE>'
  spec.setDriftTimeUnit(pyopenms.DriftTimeUnit.MILLISECOND)

  dfunit = spec.getDriftTimeUnit()
  assert dfunit == pyopenms.DriftTimeUnit.MILLISECOND
  assert pyopenms.DriftTimeUnit().getMapping()[dfunit] == "MILLISECOND"
  assert spec.getDriftTimeUnitAsString() == 'ms'

  spec = pyopenms.MSSpectrum()
  spec.setDriftTime(6.0)
  assert spec.getDriftTime() == 6.0

  spec = pyopenms.MSSpectrum()
  assert not spec.containsIMData()
  data = np.array( [5, 8, 42] ).astype(np.float32)
  f_da = [ pyopenms.FloatDataArray() ]
  f_da[0].set_data(data)
  f_da[0].setName("Ion Mobility")
  spec.setFloatDataArrays( f_da )
  assert spec.containsIMData()
  assert spec.getIMData()[0] == 0
  f_da = [ pyopenms.FloatDataArray(), pyopenms.FloatDataArray() ]
  f_da[0].setName("test")
  f_da[0].set_data(data)
  f_da[1].set_data(data)
  f_da[1].setName("Ion Mobility")
  spec.setFloatDataArrays( f_da )
  assert spec.containsIMData()
  assert spec.getIMData()[0] == 1

  # Ensure that "set_peaks()" doesnt clear the float data arrays
  spec = pyopenms.MSSpectrum()
  data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
  data_i = np.array( [50.0, 80.0] ).astype(np.float32)
  f_da = [ pyopenms.FloatDataArray() ]
  f_da[0].set_data(data)
  f_da[0].setName("Ion Mobility")
  spec.setFloatDataArrays( f_da )
  spec.set_peaks( [data_mz,data_i] )
  assert spec.containsIMData()
  assert spec.getIMData()[0] == 0
  assert len(spec.getFloatDataArrays()) == 1

  f = spec.getFloatDataArrays()[0]
  assert len(f.get_data()) == 3
  assert f.get_data()[0] == 5
  assert spec.size() == len(data_mz)
  assert spec.size() == len(data_i)

  # Test __repr__ and __str__ methods
  spec_repr = pyopenms.MSSpectrum()
  spec_repr.setRT(1234.5)
  spec_repr.setMSLevel(2)
  spec_repr.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 500.0]))

  repr_str = repr(spec_repr)
  assert "MSSpectrum(" in repr_str
  assert "ms_level=" in repr_str
  assert "rt=" in repr_str
  assert "num_peaks=" in repr_str
  str_str = str(spec_repr)
  assert str_str == repr_str


@report
def testMSChromatogram():
  """
  @tests: MSChromatogram
   MSChromatogram.__init__
   MSChromatogram.__copy__
   """
  chrom = pyopenms.MSChromatogram()
  chrom_ = copy.copy(chrom)
  assert chrom_ == chrom
  chrom_ = copy.deepcopy(chrom)
  assert chrom_ == chrom
  chrom_ = pyopenms.MSChromatogram(chrom)
  assert chrom_ == chrom

  _testMetaInfoInterface(chrom)

  chrom.setName("chrom")
  assert chrom.getName() == "chrom"

  p = pyopenms.ChromatogramPeak()
  p.setRT(1000.0)
  p.setIntensity(200.0)

  chrom.push_back(p)
  assert chrom.size() == 1
  assert chrom[0] == p

  chrom.updateRanges()
  assert isinstance(chrom.findNearest(0.0), int)

  assert isinstance(chrom.getMinRT(), float)
  assert isinstance(chrom.getMaxRT(), float)
  assert isinstance(chrom.getMinIntensity(), float)
  assert isinstance(chrom.getMaxIntensity(), float)

  assert chrom == chrom
  assert not chrom != chrom

  mz, ii = chrom.get_peaks()
  assert len(mz) == len(ii)
  assert len(mz) == 1

  chrom.set_peaks((mz, ii))
  mz0, ii0 = chrom.get_peaks()
  assert mz0 == mz
  assert ii0 == ii

  assert int(chrom.isSorted()) in (0,1)

  chrom.clear(False)
  p = pyopenms.ChromatogramPeak()
  p.setRT(1000.0)
  p.setIntensity(200.0)
  chrom.push_back(p)
  p = pyopenms.ChromatogramPeak()
  p.setRT(2000.0)
  p.setIntensity(400.0)
  chrom.push_back(p)

  mz, ii = chrom.get_peaks()
  assert chrom[0].getRT() == 1000.0
  assert chrom[1].getRT() == 2000.0
  assert chrom[0].getIntensity() == 200.0
  assert chrom[1].getIntensity() == 400.0
  assert mz[0] == 1000.0
  assert mz[1] == 2000.0
  assert ii[0] == 200.0
  assert ii[1] == 400.0

  chrom.setNativeID('chrom_0')
  chrom.setMetaValue('FWHM', 5.0)
  prec = pyopenms.Precursor()
  prec.setMZ(100.0)
  prec.setCharge(1)
  chrom.setPrecursor(prec)
  prod = pyopenms.Product()
  prod.setMZ(50.0)
  chrom.setProduct(prod)
  chrom.setComment('comment')

  df = chrom.get_df()
  # Default columns: rt, intensity, precursor_mz, precursor_charge, product_mz, native_id, FWHM
  # chromatogram_type and comment are NOT included by default
  assert df.shape == (2, 7)
  assert df.loc[0, 'rt'] == 1000.0
  assert df.loc[1, 'intensity'] == 400
  assert df.loc[1, 'precursor_mz'] == 100.0
  assert df.loc[0, 'precursor_charge'] == 1
  assert df.loc[1, 'product_mz'] == 50.0
  assert df.loc[1, 'native_id'] == 'chrom_0'
  assert df.loc[0, 'FWHM'] == 5

  # Test non-default columns (chromatogram_type, comment) via explicit selection
  df_all = chrom.get_df(columns=['rt', 'intensity', 'chromatogram_type', 'comment'])
  assert df_all.shape == (2, 4)
  assert df_all.loc[0, 'chromatogram_type'] == 'MASS_CHROMATOGRAM'
  assert df_all.loc[0, 'comment'] == 'comment'

  chrom.clear(False)
  data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
  data_i = np.array( [50.0, 80.0] ).astype(np.float32)
  chrom.set_peaks( [data_mz,data_i] )

  mz, ii = chrom.get_peaks()
  assert chrom[0].getRT() == 5.0
  assert chrom[1].getRT() == 8.0
  assert chrom[0].getIntensity() == 50.0
  assert chrom[1].getIntensity() == 80.0
  assert mz[0] == 5.0
  assert mz[1] == 8.0
  assert ii[0] == 50.0
  assert ii[1] == 80.0

  # Fast
  chrom.clear(False)
  data_mz = np.array( [5.0, 8.0] ).astype(np.float64)
  data_i = np.array( [50.0, 80.0] ).astype(np.float64)
  chrom.set_peaks( [data_mz,data_i] )

  mz, ii = chrom.get_peaks()
  assert chrom[0].getRT() == 5.0
  assert chrom[1].getRT() == 8.0
  assert chrom[0].getIntensity() == 50.0
  assert chrom[1].getIntensity() == 80.0
  assert mz[0] == 5.0
  assert mz[1] == 8.0
  assert ii[0] == 50.0
  assert ii[1] == 80.0

  # Slow
  chrom.clear(False)
  data_mz = np.array( [5.0, 8.0] ).astype(np.float32)
  data_i = np.array( [50.0, 80.0] ).astype(np.float32)
  chrom.set_peaks( [data_mz,data_i] )

  mz, ii = chrom.get_peaks()
  assert chrom[0].getRT() == 5.0
  assert chrom[1].getRT() == 8.0
  assert chrom[0].getIntensity() == 50.0
  assert chrom[1].getIntensity() == 80.0
  assert mz[0] == 5.0
  assert mz[1] == 8.0
  assert ii[0] == 50.0
  assert ii[1] == 80.0

  # test float data
  chrom = pyopenms.MSChromatogram()
  data = np.array( [5, 8, 42] ).astype(np.float32)
  f_da = [ pyopenms.FloatDataArray() ]
  f_da[0].set_data(data)
  f_da[0].setName("Test Data")
  chrom.setFloatDataArrays( f_da )
  assert len(chrom.getFloatDataArrays()) == 1
  f = chrom.getFloatDataArrays()[0]
  assert f.get_data()[0] == 5
  assert f.getName() == "Test Data"

  # Ensure that "set_peaks()" doesnt clear the float data arrays
  chrom = pyopenms.MSChromatogram()
  chrom.setFloatDataArrays( f_da )
  chrom.set_peaks( [data_mz,data_i] )
  assert len(chrom.getFloatDataArrays()) == 1

  f = chrom.getFloatDataArrays()[0]
  assert len(f.get_data()) == 3
  assert f.get_data()[0] == 5
  assert chrom.size() == len(data_mz)
  assert chrom.size() == len(data_i)

  # Test __repr__ and __str__ methods
  chrom_repr = pyopenms.MSChromatogram()
  chrom_repr.setName("test_chromatogram")
  chrom_repr.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 500.0]))

  repr_str = repr(chrom_repr)
  assert "MSChromatogram(" in repr_str
  assert "num_peaks=" in repr_str
  str_str = str(chrom_repr)
  assert str_str == repr_str


@report
def testBase64():
  """
  """

  b = pyopenms.Base64()
  out = pyopenms.String()
  inp = [1.0, 2.0, 3.0]
  b.encode64(inp, b.ByteOrder.BYTEORDER_LITTLEENDIAN, out, False)
  res = out.toString()
  assert len(res) != 0
  assert res != ""

  convBack = []
  b.decode64(res, b.ByteOrder.BYTEORDER_LITTLEENDIAN, convBack, False)
  assert convBack == inp, convBack

  # For 32 bit
  out = pyopenms.String()
  b.encode32(inp, b.ByteOrder.BYTEORDER_LITTLEENDIAN, out, False)
  res = out.toString()
  assert len(res) != 0
  assert res != ""

  convBack = []
  b.decode32(res, b.ByteOrder.BYTEORDER_LITTLEENDIAN, convBack, False)
  assert convBack == inp, convBack

@report
def testNumpressCoder():
  """
  """

  np = pyopenms.MSNumpressCoder()

  nc = pyopenms.NumpressConfig()
  nc.np_compression = np.NumpressCompression.LINEAR
  nc.estimate_fixed_point = True
  tmp = pyopenms.String()
  out = []
  inp = [1.0, 2.0, 3.0]
  np.encodeNP(inp, tmp, True, nc)

  res = tmp.toString()
  assert len(res) != 0, len(res)
  assert res != "", res
  np.decodeNP(res, out, True, nc)
  assert len(out) == 3, (out, res)
  assert out == inp, out

  # Now try to use a simple Python string as input -> this will fail as we
  # cannot pass this by reference in C++
  res = ""
  try:
    np.encodeNP(inp, res, True, nc)
    has_error = False
  except AssertionError:
    has_error = True

  assert has_error

@report
def testNumpressConfig():
  """
  """

  n = pyopenms.MSNumpressCoder()
  np = pyopenms.NumpressConfig()
  np.np_compression = n.NumpressCompression.LINEAR
  assert np.np_compression == n.NumpressCompression.LINEAR
  np.numpressFixedPoint = 4.2
  np.numpressErrorTolerance = 4.2
  np.estimate_fixed_point = True
  np.linear_fp_mass_acc = 4.2
  np.setCompression("linear")

@report
def testPeakTypeEstimator():
  """
  @tests: PeakTypeEstimator
   PeakTypeEstimator.__init__
   PeakTypeEstimator.estimateType
  """

  pyopenms.PeakTypeEstimator().estimateType(pyopenms.MSSpectrum())


if __name__ == "__main__":
  # Run all tests in this module
  import sys
  module = sys.modules[__name__]
    
  test_functions = [getattr(module, name) for name in dir(module) 
           if name.startswith('test') and callable(getattr(module, name))]
    
  print(f"Running {len(test_functions)} spectrum/experiment tests...")
  for test_func in test_functions:
    try:
      test_func()
      print(f"✅ {test_func.__name__}")
    except (AssertionError, RuntimeError, ValueError, TypeError) as e:
      print(f"❌ {test_func.__name__}: {e}")
