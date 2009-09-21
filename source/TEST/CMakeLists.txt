### Build rules for Class tests
PROJECT("OpenMS_tests")

enable_testing()
INCLUDE(Dart)
	
## Warning: this needs to be the same as the TOPP tools dir (where the OpenMS.dll will be on windows) - otherwise the tests won't run!
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
include(executables.cmake)

if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS_RELEASE_SAVED ${CMAKE_CXX_FLAGS_RELEASE})
	#disable optimization for tests
	set(CMAKE_CXX_FLAGS_RELEASE "-O0")
endif()


### CONCEPT
add_custom_target(CONCEPT_TEST)
add_dependencies(CONCEPT_TEST ${concept_executables_list})

### DATASTRUCTURES
add_custom_target(DATASTRUCTURES_TEST)
add_dependencies(DATASTRUCTURES_TEST ${datastructures_executables_list})

### METADATA
add_custom_target(METADATA_TEST)
add_dependencies(METADATA_TEST ${metadata_executables_list})

### SYSTEM
add_custom_target(SYSTEM_TEST)
add_dependencies(SYSTEM_TEST ${system_executables_list})

### KERNEL
add_custom_target(KERNEL_TEST)
add_dependencies(KERNEL_TEST ${kernel_executables_list})

### VISUAL
add_custom_target(VISUAL_TEST)
add_dependencies(VISUAL_TEST ${visual_executables_list})

### FORMAT
add_custom_target(FORMAT_TEST)
add_dependencies(FORMAT_TEST ${format_executables_list})

### MATH
add_custom_target(MATH_TEST)
add_dependencies(MATH_TEST ${math_executables_list})

### FILTERING
add_custom_target(FILTERING_TEST)
add_dependencies(FILTERING_TEST ${filtering_executables_list})

### COMPARISON
add_custom_target(COMPARISON_TEST)
add_dependencies(COMPARISON_TEST ${comparison_executables_list})

### CHEMISTRY
add_custom_target(CHEMISTRY_TEST)
add_dependencies(CHEMISTRY_TEST ${chemistry_executables_list})

### ANALYSIS
add_custom_target(ANALYSIS_TEST)
add_dependencies(ANALYSIS_TEST ${analysis_executables_list})

### APPLICATIONS
add_custom_target(APPLICATIONS_TEST)
add_dependencies(APPLICATIONS_TEST ${applications_executables_list})

### TRANSFORMATIONS
add_custom_target(TRANSFORMATIONS_TEST)
add_dependencies(TRANSFORMATIONS_TEST ${transformations_executables_list})

### SIMULATION
add_custom_target(SIMULATION_TEST)
add_dependencies(SIMULATION_TEST ${simulation_executables_list})

if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX)
	set(${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_RELEASE_SAVED})
endif()

										
### tests
foreach(i ${TEST_executables})
	add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
	# should not be necessary if ctest_build entries are also done for this TEST directory
	#add_dependencies(${i} TOPP UTILS Tutorials_build FuzzyDiff)
	add_test(${i} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
endforeach(i)

### "test_build" target
add_custom_target(test_build)
# UTILS and Tutorials_build are not necessarily needed, however very useful
# for nightly builds (untested otherwise!)
add_dependencies(test_build ${TEST_executables} TOPP UTILS Tutorials_build FuzzyDiff)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TEST_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})

############## TOPP test ####################
add_subdirectory(TOPP)
