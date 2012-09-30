project(OpenSwath)

set(OpenSwathAlgorithmFiles
  ALGO/Scoring.C
  ALGO/MRMScoring.C
  ALGO/DIAHelpers.C
)

# message(STATUS ">>>>>>openswath<<<<<<")
# message(STATUS "src -------- ${src} -------")
# message(STATUS "lib -------- ${lib} -------")
# message(STATUS "inc -------- ${inc} -------")

include_directories(${inc})
##link_directories(${lib})
## library with additional classes from above
add_library(OpenSwathAlgo SHARED ${OpenSwathAlgorithmFiles})

