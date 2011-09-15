### Build rules for Class tests
PROJECT("OpenMS_tests") ## do not change this name, as you will break Windows-NightlyBuilds

enable_testing()
INCLUDE(Dart)

# required to build GUI tests (e.g., TOPPView_test)
ADD_DEFINITIONS(-DQT_GUI_LIB)


## - BEGIN OF WORKAROUND
## this is a workaround for a CMake-Bug (http://www.cmake.org/Bug/view.php?id=12419)
## , which prevents the orginal svn_revision_update target to be build prior to OpenMS.lib on Win-VS
## Delete when bug is fixed.
##
add_custom_target(svn_revision_update_forTests ALL
	    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -DSVN_REVISION_FILE=${SVN_REVISION_FILE} -DOPENMS_HAS_SVNVERSION=${OPENMS_HAS_SVNVERSION} -DSVNVERSION_EXECUTABLE=${SVNVERSION_EXECUTABLE}
	                             -P ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/svnversion.cmake)
add_dependencies(OpenMS svn_revision_update_forTests)
## - END OF WORKAROUND

## Warning: in this directory the OpenMS.dll should already be present (e.g. by using the Dll's POST_BUILD command to copy it)  
##          - otherwise the tests might not run!
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
include(executables.cmake)

if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
	#disable optimization for tests
	set(CMAKE_CXX_FLAGS_RELEASE_SAVED ${CMAKE_CXX_FLAGS_RELEASE})
	set(CMAKE_CXX_FLAGS_RELEASE "-O0")
	if(MT_ENABLE_OPENMP AND ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	  # disable OPENMP for Mac OS X/gcc, see
	  # https://sourceforge.net/apps/trac/open-ms/ticket/280
	  # bug for details
	  STRING(REPLACE "-fopenmp" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
	endif()
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

### VISUAL
add_custom_target(VISUAL_TEST)
add_dependencies(VISUAL_TEST ${visual_executables_list})

### tests
foreach(i ${TEST_executables})
	add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
	add_test(${i} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
	# only add OPENMP flags to gcc linker (execpt Mac OS X, due to compiler bug
	# see https://sourceforge.net/apps/trac/open-ms/ticket/280 for details)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()  
endforeach(i)

### tests of classes from GUI lib:
foreach(i ${visual_executables_list})
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
endforeach(i)

## GUI tests special treatment
### Apply MOC compiler to GUI
set(GUI_executables_list_C)
foreach(i ${GUI_executables_list})
  list(APPEND GUI_executables_list_C GUI/${i}.C)
endforeach(i)
## create mocced files - they are included manually in the implementation files
QT4_WRAP_CPP(mocced_GUI_executables_list_C ${GUI_executables_list_C})

foreach(i ${GUI_executables_list})
	SET_SOURCE_FILES_PROPERTIES(GUI/${i}.C PROPERTIES OBJECT_DEPENDS ${mocced_GUI_executables_list_C})
	add_executable(${i} GUI/${i}.C)
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	add_test(${i} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
	# only add OPENMP flags to gcc linker (execpt Mac OS X, due to compiler bug
	# see https://sourceforge.net/apps/trac/open-ms/ticket/280 for details)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()  
endforeach(i)

### GUI tests
add_custom_target(GUI_TEST)
add_dependencies(GUI_TEST ${GUI_executables_list})

### add gui tests to executables list
set(TEST_executables ${TEST_executables} TOPPView_test)

if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
	set(${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_RELEASE_SAVED})
endif()

### "test_build" target
add_custom_target(test_build ALL)
add_dependencies(test_build ${TEST_executables} TOPP UTILS)

### all target that we require to show up in the targets list
### and in the nightly build log (to test if they can be build properly)
add_custom_target(extra_targets ALL)
add_dependencies(extra_targets Tutorials_build)


### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TEST_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})

############## TOPP test ####################
add_subdirectory(TOPP)
