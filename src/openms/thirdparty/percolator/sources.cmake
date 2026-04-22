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
)
list(TRANSFORM OpenMS_percolator_sources PREPEND "${CMAKE_CURRENT_LIST_DIR}/")
