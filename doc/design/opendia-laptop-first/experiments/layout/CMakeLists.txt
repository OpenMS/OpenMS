cmake_minimum_required(VERSION 3.24)
project(OpenDIALayoutExperiment LANGUAGES CXX)
# $Maintainer: Justin Sing $
# Build out of tree. The WP01 checkout is read-only input, not a vendored copy.
set(OPENDIA_SOURCE_DIR "" CACHE PATH "Required read-only WP01 checkout at 3cd24a2502")
if(NOT EXISTS "${OPENDIA_SOURCE_DIR}/src/topp/OpenDIA/source/OpenDIARecords.cpp")
  message(FATAL_ERROR "Set OPENDIA_SOURCE_DIR to the WP01 source checkout")
endif()
find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
add_executable(layout_benchmark layout_benchmark.cpp
  "${OPENDIA_SOURCE_DIR}/src/topp/OpenDIA/source/OpenDIARecords.cpp"
  "${OPENDIA_SOURCE_DIR}/src/openms/source/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.cpp")
target_include_directories(layout_benchmark PRIVATE "${OPENDIA_SOURCE_DIR}/src/topp/OpenDIA/include")
target_link_libraries(layout_benchmark PRIVATE OpenMS::OpenMS_CLI)
target_compile_features(layout_benchmark PRIVATE cxx_std_23)
# Apply identical optimization to benchmark and actual production record implementation.
target_compile_options(layout_benchmark PRIVATE -O3 -DNDEBUG -Wall -Wextra -Wpedantic)
include(CTest)
if(BUILD_TESTING)
  add_test(NAME layout_self_test COMMAND layout_benchmark --self-test)
  foreach(variant rawtile aos soa factored-aos factored-soa)
    foreach(fixture topp/OpenSwathWorkflow_1_input topp/OpenSwathWorkflow_23_input class_tests/openms/data/IM_FAIMS_test class_tests/openms/data/ChromatogramExtractor_input)
      string(REPLACE "/" "_" test_name "${fixture}")
      add_test(NAME "layout_${variant}_${test_name}" COMMAND layout_benchmark
        "${OPENDIA_SOURCE_DIR}/src/tests/${fixture}.mzML" "${variant}" 0 1)
    endforeach()
  endforeach()
endif()
