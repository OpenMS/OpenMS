# src/openms/thirdparty/percolator/sources.cmake
# Lists the .cpp files to be compiled as part of libOpenMS.

set(OpenMS_percolator_sources
  BaseSpline.cpp
  CrossValidation.cpp
  DataSet.cpp
  FeatureMemoryPool.cpp
  FeatureNames.cpp
  Globals.cpp
  IsotonicPEP.cpp
  Logger.cpp
  LogisticRegression.cpp
  MassHandler.cpp
  MonotoneRegressor.cpp
  MyException.cpp
  NoNormalizer.cpp
  Normalizer.cpp
  Numerical.cpp
  PackedMatrix.cpp
  PackedVector.cpp
  PosteriorEstimator.cpp
  PSMDescription.cpp
  PseudoRandom.cpp
  ResultHolder.cpp
  SanityCheck.cpp
  ScoreHolder.cpp
  Scores.cpp
  Set.cpp
  SetHandler.cpp
  StdvNormalizer.cpp
  Timer.cpp
  UniNormalizer.cpp
  Vector.cpp
  ssl.cpp

  # OpenMS-local addition: reservoir-sampling helper derived from
  # SetHandler::readPSMs, required by the in-process Percolator wrapper
  # since we don't vendor Caller.cpp. Not in whitelist.txt — preserved
  # across sync-from-upstream.sh runs. See InProcessSampler.h for details.
  InProcessSampler.cpp
)
list(TRANSFORM OpenMS_percolator_sources PREPEND "${CMAKE_CURRENT_LIST_DIR}/")
