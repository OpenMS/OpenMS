cmake_minimum_required(VERSION 3.24)
project(native_cid LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
get_filename_component(ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(BUNDLE "${ROOT}/../prose-andes-reproduction")
set(OPENMS "${ROOT}/sources/base" CACHE PATH "OpenMS source checkout")
set(TESTCFG "" CACHE PATH "Directory with a per-build OpenMS/test_config.h")
include_directories(BEFORE "${TESTCFG}")
# ONNX variant: PeptDeep rescoring (#9975) compiled with ONNX Runtime 1.23.2 (develop's vcpkg pin)
set(ORT "${ROOT}/../onnx/onnxruntime-linux-x64-1.23.2")
add_compile_definitions(WITH_ONNX=1)
include_directories("${ORT}/include" "${OPENMS}/src/openms/source/ML/ONNX")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include" "${OPENMS}/src/openms/include" "${OPENMS}/src/openswathalgo/include" "${OPENMS}/src/openms/extern/nlohmann_json/include" "${BUNDLE}/header_deps/cmeel.prefix/include" "${BUNDLE}/header_deps/cmeel.prefix/include/eigen3")
add_library(prose_objects OBJECT
 "${OPENMS}/src/openms/source/ANALYSIS/ID/ProSEAlgorithm.cpp"
 "${OPENMS}/src/openms/source/ANALYSIS/ID/FragmentIndex.cpp"
 "${OPENMS}/src/openms/source/ANALYSIS/ID/HyperScore.cpp"
 "${OPENMS}/src/openms/source/CHEMISTRY/TheoreticalSpectrumGenerator.cpp"
 "${OPENMS}/src/openms/source/ANALYSIS/ID/PeptDeepRescoring.cpp"
 "${OPENMS}/src/openms/source/ML/ONNX/ONNXEnvironment.cpp"
 "${OPENMS}/src/openms/source/ML/ONNX/ONNXPredictorBase.cpp"
 "${OPENMS}/src/openms/source/ML/PEPTDEEP/PeptDeepInput.cpp"
 "${OPENMS}/src/openms/source/ML/PEPTDEEP/PeptDeepMS2Inference.cpp"
 "${OPENMS}/src/openms/source/ML/PEPTDEEP/PeptDeepRTInference.cpp"
 "${OPENMS}/src/openms/source/ML/PEPTDEEP/PeptDeepCCSInference.cpp")
if(EXISTS "${OPENMS}/src/openms/source/ANALYSIS/ID/FragmentIonLikelihoodModel.cpp")
  target_sources(prose_objects PRIVATE "${OPENMS}/src/openms/source/ANALYSIS/ID/FragmentIonLikelihoodModel.cpp")
endif()
target_compile_options(prose_objects PRIVATE -fopenmp)
file(GLOB OPENMS_LIB "${BUNDLE}/nightly/pyopenms.libs/libOpenMS-*.so")
add_executable(prose main.cpp $<TARGET_OBJECTS:prose_objects>)
target_link_libraries(prose PRIVATE "${OPENMS_LIB}" "${ORT}/lib/libonnxruntime.so" -fopenmp)
target_link_options(prose PRIVATE "-Wl,-rpath-link,${BUNDLE}/nightly/pyopenms.libs")
set_target_properties(prose PROPERTIES BUILD_RPATH "${BUNDLE}/nightly/pyopenms.libs;${ORT}/lib")
include_directories("${OPENMS}/src/testframework/include")
add_library(test_support OBJECT
  "${OPENMS}/src/testframework/source/CONCEPT/ClassTest.cpp"
  "${OPENMS}/src/testframework/source/CONCEPT/FuzzyStringComparator.cpp"
  "${OPENMS}/src/tests/class_tests/source/OpenMSTestSupport.cpp")
enable_testing()
set(TESTS HyperScore ProSEAlgorithm FragmentIndex PeptDeepRescoring PeptDeepInference)
if(EXISTS "${OPENMS}/src/tests/class_tests/openms/source/FragmentIonLikelihoodModel_test.cpp")
  list(APPEND TESTS FragmentIonLikelihoodModel)
endif()
foreach(TEST ${TESTS})
  add_executable(${TEST}_test "${OPENMS}/src/tests/class_tests/openms/source/${TEST}_test.cpp" $<TARGET_OBJECTS:prose_objects> $<TARGET_OBJECTS:test_support>)
  target_link_libraries(${TEST}_test PRIVATE "${OPENMS_LIB}" "${ORT}/lib/libonnxruntime.so" -fopenmp)
  target_link_options(${TEST}_test PRIVATE "-Wl,-rpath-link,${BUNDLE}/nightly/pyopenms.libs")
  set_target_properties(${TEST}_test PROPERTIES BUILD_RPATH "${BUNDLE}/nightly/pyopenms.libs;${ORT}/lib")
  add_test(NAME ${TEST}_test COMMAND ${TEST}_test WORKING_DIRECTORY "${ROOT}/../onnx/testrun")
endforeach()
