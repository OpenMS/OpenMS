project(tests)

set(my_tests
mQuestScoring_test
Scoring_test
TestConvert
DiaHelpers_test
)

## add targets for the executables
foreach(i ${my_tests})
  add_executable(${i} ${i}.C)
  target_link_libraries(${i} OpenSwathAlgo boost_unit_test_framework)
	add_test(${i} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
endforeach(i)

