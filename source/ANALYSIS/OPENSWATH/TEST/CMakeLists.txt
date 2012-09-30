project(tests)

set(my_tests
mQuestScoring_test
Scoring_test
TestConvert
DiaHelpers_test
)

message( STATUS ">>>>>>teststests<<<<<<" )
message( STATUS "src -------- ${src} -------" )
message( STATUS "lib -------- ${lib} -------" )
message( STATUS "inc -------- ${inc} -------" )

include_directories(${inc})
link_directories(${lib})

## add targets for the executables
foreach(i ${my_tests})
  add_executable(${i} ${i}.C)
  target_link_libraries(${i} OpenSwathAlgo boost_unit_test_framework)
endforeach(i)

ADD_TEST(NAME mQuestScoring_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} COMMAND COMMAND ${BinOut}/mQuestScoring_test)
ADD_TEST(NAME Scoring_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} COMMAND COMMAND ${BinOut}/Scoring_test)
ADD_TEST(NAME DiaHelpers_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} COMMAND COMMAND ${BinOut}/DiaHelpers_test)
ENABLE_TESTING()
