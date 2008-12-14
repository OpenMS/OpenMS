##

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/source/TEST)
include(executables.cmake)

### CONCEPT
foreach(i ${concept_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(CONCEPT_TEST)
add_dependencies(CONCEPT_TEST ${concept_executables_list})

### DATASTRUCTURES
foreach(i ${datastructures_executables_list})
	add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(DATASTRUCTURES_TEST)
add_dependencies(DATASTRUCTURES_TEST ${datastructures_executables_list})

### METADATA
foreach(i ${metadata_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(METADATA_TEST)
add_dependencies(METADATA_TEST ${metadata_executables_list})

### SYSTEM
foreach(i ${system_executables_list})
  add_executable(${i} ${i}.C)
  target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(SYSTEM_TEST)
add_dependencies(SYSTEM_TEST ${system_executables_list})

### KERNEL
foreach(i ${kernel_executables_list})
  add_executable(${i} ${i}.C)
  target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(KERNEL_TEST)
add_dependencies(KERNEL_TEST ${kernel_executables_list})


### VISUAL
foreach(i ${visual_executables_list})
  add_executable(${i} ${i}.C)
  target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(VISUAL_TEST)
add_dependencies(VISUAL_TEST ${visual_executables_list})

### FORMAT
foreach(i ${format_executables_list})
  add_executable(${i} ${i}.C)
  target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(FORMAT_TEST)
add_dependencies(FORMAT_TEST ${format_executables_list})

### MATH
foreach(i ${math_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(MATH_TEST)
add_dependencies(MATH_TEST ${math_executables_list})

### FILTERING
foreach(i ${filtering_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(FILTERING_TEST)
add_dependencies(FILTERING_TEST ${filtering_executables_list})

### COMPARISON
foreach(i ${comparison_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(COMPARISON_TEST)
add_dependencies(COMPARISON_TEST ${comparison_executables_list})

### CHEMISTRY
foreach(i ${chemistry_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(CHEMISTRY_TEST)
add_dependencies(CHEMISTRY_TEST ${chemistry_executables_list})

### ANALYSIS
foreach(i ${analysis_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(ANALYSIS_TEST)
add_dependencies(ANALYSIS_TEST ${analysis_executables_list})
		
### APPLICATIONS
foreach(i ${applications_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(APPLICATIONS_TEST)
add_dependencies(APPLICATIONS_TEST ${applications_executables_list})

### TRANSFORMATIONS
foreach(i ${transformations_executables_list})
  add_executable(${i} ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(TRANSFORMATIONS_TEST)
add_dependencies(TRANSFORMATIONS_TEST ${transformations_executables_list})


### collect test executables
set(TEST_executables 
		${concept_executables_list} 
		${system_executables_list} 
		${datastructures_executables_list} 
		${kernel_executables_list} 
		${metadata_executables_list}
		${visual_executables_list}
		${format_executables_list}
		${math_executables_list}
		${filtering_executables_list}
		${comparison_executables_list}
		${chemistry_executables_list}
		${analysis_executables_list}
		${applications_executables_list}
		${transformations_executables_list})

### TEST target 
add_custom_target(TEST)
add_dependencies(TEST ${TEST_executables})

### tests
foreach(i ${TEST_executables})
  add_test(${i} ${i})
endforeach(i)

