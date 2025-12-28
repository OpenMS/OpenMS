#!/usr/bin/env python
# -*- coding: utf-8  -*-

## ----------------------------------------------------------------------------
## Core Algorithms tests extracted from test000.py
## Part of Issue #8567: Split test000.py into modular test files
## ----------------------------------------------------------------------------

from __future__ import print_function
import pyopenms
import numpy as np

# Import shared test helper functions
from test_helpers import report, _testParam, _testStrOutput


@report
def testMapAlignment():
  """
  @tests: MapAlignmentAlgorithmPoseClustering
   MapAlignmentAlgorithmPoseClustering.__init__
   MapAlignmentAlgorithmPoseClustering.getDefaults
   MapAlignmentAlgorithmPoseClustering.getName
   MapAlignmentAlgorithmPoseClustering.getParameters
   MapAlignmentAlgorithmPoseClustering.setName
   MapAlignmentAlgorithmPoseClustering.setParameters
   MapAlignmentAlgorithmPoseClustering.setReference
   MapAlignmentAlgorithmPoseClustering.align
  """
  ma = pyopenms.MapAlignmentAlgorithmPoseClustering()
  assert isinstance(ma.getDefaults(), pyopenms.Param)
  assert isinstance(ma.getParameters(), pyopenms.Param)
  _testStrOutput(ma.getName())

  ma.setName(ma.getName())


  assert ma.getDefaults() is not None
  assert ma.getParameters() is not None

  ma.setParameters(ma.getDefaults())

  assert ma.setReference is not None
  assert ma.align is not None


@report
def testMapAlignmentIdentification():
  """
  @tests: MapAlignmentAlgorithmIdentification
   MapAlignmentAlgorithmIdentification.__init__
   MapAlignmentAlgorithmIdentification.align
   MapAlignmentAlgorithmIdentification.setReference
   """
  ma = pyopenms.MapAlignmentAlgorithmIdentification()

  assert pyopenms.MapAlignmentAlgorithmIdentification().align is not None
  assert pyopenms.MapAlignmentAlgorithmIdentification().setReference is not None


@report
def testMapAlignmentTransformer():
  """
  @tests: MapAlignmentTransformer
   MapAlignmentTransformer.__init__
   MapAlignmentTransformer.transformRetentionTimes
   """
  ma = pyopenms.MapAlignmentTransformer()

  assert pyopenms.MapAlignmentTransformer().transformRetentionTimes is not None


@report
def testTransformationDescription():
  """
  @tests: TransformationDescription
   TransformationDescription.__init__
   TransformationDescription.apply
   TransformationDescription.getDataPoints
   TransformationDescription.fitModel
   TransformationDescription.getModelParameters
   TransformationDescription.getModelType
   TransformationDescription.invert
  """
  td = pyopenms.TransformationDescription()
  assert td.getDataPoints() == []
  assert isinstance(td.apply(0.0), float)

  td.fitModel
  p = td.getModelParameters()
  td.getModelType()
  td.invert


@report
def testTransformationModels():
  """
  @tests: TransformationModelInterpolated
   TransformationModelInterpolated.getDefaultParameters
   TransformationModelInterpolated.getParameters
   TransformationModelLinear.getDefaultParameters
   TransformationModelLinear.getParameters
   TransformationModelBSpline.getDefaultParameters
   TransformationModelBSpline.getParameters
   TransformationModelLowess.getDefaultParameters
   TransformationModelLowess.getParameters
  """
  for clz in [pyopenms.TransformationModelLinear,
        pyopenms.TransformationModelBSpline,
        pyopenms.TransformationModelInterpolated,
        pyopenms.TransformationModelLowess]:
    p = pyopenms.Param()
    data = [ pyopenms.TM_DataPoint(9.0, 8.9),
         pyopenms.TM_DataPoint(5.0, 6.0),
         pyopenms.TM_DataPoint(8.0, 8.0) ]
    mod = clz(data, p)
    mod.evaluate(7.0)
    mod.getDefaultParameters(p)


@report
def testLinearResampler():
  """
  @tests: LinearResampler
   LinearResampler.__init__
   LinearResampler.raster
   LinearResampler.rasterExperiment
  """
  ff = pyopenms.LinearResampler()
  p = ff.getDefaults()
  _testParam(p)

  assert pyopenms.LinearResampler().raster is not None
  assert pyopenms.LinearResampler().rasterExperiment is not None


@report
def testPeakPickerChromatogram():
  """
  @tests: PeakPickerChromatogram
   PeakPickerChromatogram.__init__
   PeakPickerChromatogram.pickChromatogram
  """
  p = pyopenms.PeakPickerChromatogram()
  assert p.pickChromatogram is not None


@report
def testPeakPickerHiRes():
  """
  @tests: PeakPickerHiRes
   PeakPickerHiRes.__init__
   PeakPickerHiRes.endProgress
   PeakPickerHiRes.getDefaults
   PeakPickerHiRes.getLogType
   PeakPickerHiRes.getName
   PeakPickerHiRes.getParameters
   PeakPickerHiRes.pick
   PeakPickerHiRes.pickExperiment
   PeakPickerHiRes.setLogType
   PeakPickerHiRes.setName
   PeakPickerHiRes.setParameters
   PeakPickerHiRes.setProgress
   PeakPickerHiRes.startProgress
  """
  p = pyopenms.PeakPickerHiRes()
  assert p.pick is not None
  assert p.pickExperiment is not None


@report
def testConsensusMapNormalizerAlgorithmMedian():
  """
  @tests: ConsensusMapNormalizerAlgorithmMedian
   ConsensusMapNormalizerAlgorithmMedian.__init__
   ConsensusMapNormalizerAlgorithmMedian.normalizeMaps
  """
  ff = pyopenms.ConsensusMapNormalizerAlgorithmMedian()

  assert pyopenms.ConsensusMapNormalizerAlgorithmMedian().normalizeMaps is not None


@report
def testConsensusMapNormalizerAlgorithmQuantile():
  """
  @tests: ConsensusMapNormalizerAlgorithmQuantile
   ConsensusMapNormalizerAlgorithmQuantile.__init__
   ConsensusMapNormalizerAlgorithmQuantile.normalizeMaps
  """
  ff = pyopenms.ConsensusMapNormalizerAlgorithmQuantile()

  assert pyopenms.ConsensusMapNormalizerAlgorithmQuantile().normalizeMaps is not None


@report
def testConsensusMapNormalizerAlgorithmThreshold():
  """
  @tests: ConsensusMapNormalizerAlgorithmThreshold
   ConsensusMapNormalizerAlgorithmThreshold.__init__
   ConsensusMapNormalizerAlgorithmThreshold.computeCorrelation
   ConsensusMapNormalizerAlgorithmThreshold.normalizeMaps
  """
  ff = pyopenms.ConsensusMapNormalizerAlgorithmThreshold()

  assert pyopenms.ConsensusMapNormalizerAlgorithmThreshold().computeCorrelation is not None
  assert pyopenms.ConsensusMapNormalizerAlgorithmThreshold().normalizeMaps is not None


# NOTE: testBSpline2d may hang/crash in some pyOpenMS versions
# This is also the case in the original test000.py
@report
def testBSpline2d():
  """
  @tests: BSpline2d
   BSpline2d.__init__
   BSpline2d.ok
   BSpline2d.eval
   BSpline2d.derivative
   BSpline2d.solve
  """
  x = [1.0, 6.0, 8.0, 10.0, 15.0]
  y = [2.0, 5.0, 6.0, 12.0, 13.0]
  spline = pyopenms.BSpline2d(x, y, 0, pyopenms.BoundaryCondition.BC_ZERO_ENDPOINTS, 0)

  assert spline.ok()
  assert (abs(spline.eval(6.0) - 5.0) < 0.01)
  assert (abs(spline.derivative(6.0) - 5.0) < 0.01)

  y_new = [4.0, 5.0, 6.0, 12.0, 13.0]
  spline.solve(y_new)

  assert spline.ok()
  assert (abs(spline.eval(6.0) - 5.0) < 0.01)


@report
def testMatrixDouble():
  """
  @tests: MatrixDouble
   MatrixDouble.__init__
   MatrixDouble.rows
   MatrixDouble.cols
   MatrixDouble.getValue
   MatrixDouble.setValue
   MatrixDouble.get_matrix
   MatrixDouble.get_matrix_as_view
   MatrixDouble.set_matrix
   """
  m = pyopenms.MatrixDouble(3, 2, 0.0)
  for i in range(3):
    for j in range(2):
      m.setValue(i, j, i * 10.0 + j) 
  print(m)

  mv = m.get_matrix_as_view()
  print(mv)

  mc = m.get_matrix()
  print(mc)

  mat = m.get_matrix_as_view()

  N = 90
  m = pyopenms.MatrixDouble(N-1, N+2, 5.0)

  assert m.rows() == 89
  assert m.cols() == 92

  rows = N-1
  cols = N+2
  test = []
  for i in range(int(rows)):
    for j in range(int(cols)):
      test.append( m.getValue(i,j) )

  testm = np.asarray(test)
  testm = testm.reshape(rows, cols)

  assert sum(sum(testm)) == 40940.0
  assert sum(sum(testm)) == (N-1)*(N+2)*5

  matrix = m.get_matrix()
  assert sum(sum(matrix)) == 40940.0
  assert sum(sum(matrix)) == (N-1)*(N+2)*5

  matrix_view = m.get_matrix_as_view()
  assert sum(sum(matrix_view)) == 40940.0
  assert sum(sum(matrix_view)) == (N-1)*(N+2)*5


  # Column = 1 / Row = 2
  ## Now change a value:

  assert m.getValue(1, 2) == 5.0
  m.setValue(1, 2, 8.0)
  assert m.getValue(1, 2) == 8.0

  print(m)
  mat = m.get_matrix_as_view()
  print(mat)
  assert mat[1, 2] == 8.0

  mat = m.get_matrix()
  assert m.getValue(1, 2) == 8.0
  assert mat[1, 2] == 8.0

  # Whatever we change here gets changed in the raw data as well
  matrix_view = m.get_matrix_as_view()
  matrix_view[1, 6] = 11.0
  assert m.getValue(1, 6) == 11.0
  assert matrix_view[1, 6] == 11.0

  m = pyopenms.MatrixDouble()
  assert m.rows() == 0
  assert m.cols() == 0

  mat[3, 6] = 9.0
  m.set_matrix(mat)
  assert m.getValue(1, 2) == 8.0
  assert m.getValue(3, 6) == 9.0

  # Create a test matrix with known values
  test_matrix = pyopenms.MatrixDouble(3, 4, 0.0)
  test_matrix.setValue(2, 3, 120.7)

  # Test getValue retrieves correct values
  assert test_matrix.getValue(2, 3) == 120.7


@report
def testMatrixDoubleColumnMajorOrdering():
  """
  @tests: MatrixDouble
  Verify Matrix preserves row/column ordering through Python<->C++ round-trip.
  This test catches column-major vs row-major issues that could cause data
  transposition when passing matrices between numpy and C++.
  """
  # Test 1: Non-square matrix with unique values at each position
  # Using 3x4 matrix where value[i,j] = i*10 + j makes each element unique
  # and any transposition immediately detectable
  original = np.array([[0, 1, 2, 3],
             [10, 11, 12, 13],
             [20, 21, 22, 23]], dtype=np.float64)  # 3 rows x 4 cols

  m = pyopenms.MatrixDouble()
  m.set_matrix(original)

  # Verify shape preserved (not transposed)
  assert m.rows() == 3, f"Expected 3 rows, got {m.rows()}"
  assert m.cols() == 4, f"Expected 4 cols, got {m.cols()}"

  # Verify each element via getValue matches numpy indexing
  for i in range(3):
    for j in range(4):
      expected = original[i, j]
      actual = m.getValue(i, j)
      assert actual == expected, \
        f"getValue mismatch at ({i},{j}): C++={actual}, numpy={expected}"

  # Test 2: Round-trip preservation (numpy -> C++ -> numpy)
  result = m.get_matrix()
  assert result.shape == original.shape, \
    f"Shape mismatch after round-trip: {result.shape} vs {original.shape}"
  assert np.array_equal(original, result), \
    f"Data mismatch after round-trip:\nOriginal:\n{original}\nResult:\n{result}"

  # Test 3: Verify view indexing matches getValue for all elements
  view = m.get_matrix_as_view()
  assert view.shape == (3, 4), f"View shape mismatch: {view.shape}"
  for i in range(3):
    for j in range(4):
      assert view[i, j] == m.getValue(i, j), \
        f"View mismatch at ({i},{j}): view={view[i,j]}, getValue={m.getValue(i,j)}"

  # Test 4: Verify modifications through view are reflected in C++ object
  view[2, 3] = 99.0
  assert m.getValue(2, 3) == 99.0, "View modification not reflected in C++ object"

  # Test 5: Test with transposed-like access pattern to catch subtle bugs
  # Create matrix where row index and col index have very different values
  m2 = pyopenms.MatrixDouble(2, 5, 0.0)  # 2 rows, 5 cols - very non-square
  for i in range(2):
    for j in range(5):
      m2.setValue(i, j, i * 100 + j)

  view2 = m2.get_matrix_as_view()
  assert view2.shape == (2, 5), f"Non-square view shape wrong: {view2.shape}"

  # Check corners and middle to ensure no transposition
  assert view2[0, 0] == 0, f"[0,0] = {view2[0,0]}, expected 0"
  assert view2[0, 4] == 4, f"[0,4] = {view2[0,4]}, expected 4"
  assert view2[1, 0] == 100, f"[1,0] = {view2[1,0]}, expected 100"
  assert view2[1, 4] == 104, f"[1,4] = {view2[1,4]}, expected 104"


@report
def testMapConversion():
  """
  @tests: MapConversion
   MapConversion.convert
  """
  feature = pyopenms.Feature()
  feature.setRT(99)

  cmap = pyopenms.ConsensusMap()
  fmap = pyopenms.FeatureMap()
  fmap.push_back(feature)
  pyopenms.MapConversion().convert(0, fmap, cmap, 1)

  assert(cmap.size() == 1)
  assert(cmap[0].getRT() == 99.0)

  fmap = pyopenms.FeatureMap()
  pyopenms.MapConversion().convert(cmap, True, fmap)

  assert(fmap.size() == 1)


if __name__ == "__main__":
  # Run all tests in this module
  import sys
  module = sys.modules[__name__]
    
  test_functions = [getattr(module, name) for name in dir(module) 
           if name.startswith('test') and callable(getattr(module, name))]
    
  print(f"Running {len(test_functions)} algorithm tests...")
  for test_func in test_functions:
    try:
      test_func()
      print(f"✅ {test_func.__name__}")
    except Exception as e:
      print(f"❌ {test_func.__name__}: {e}")
