## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "OpenMS")
set(CTEST_NIGHTLY_START_TIME "23:59:59 UTC")
set(CTEST_SUBMIT_URL "https://cdash.seqan.de/submit.php?project=OpenMS")
