cmake_minimum_required(VERSION 3.16)
project(DateTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Find Qt5
find_package(Qt5Core REQUIRED)

# Find Boost
find_package(Boost 1.74.0 REQUIRED)

# Define OpenMS version info
add_definitions(-DOPENMS_VERSION="3.6.0")

# Include directories
include_directories(
  ${CMAKE_SOURCE_DIR}/../src/openms/include
)

# Create config header
configure_file(
  ${CMAKE_SOURCE_DIR}/../src/openms/include/OpenMS/OpenMSConfig.h.in
  ${CMAKE_BINARY_DIR}/include/OpenMS/OpenMSConfig.h
)
include_directories(${CMAKE_BINARY_DIR}/include)

# Add test config
configure_file(
  ${CMAKE_SOURCE_DIR}/../src/tests/class_tests/openms/include/OpenMS/test_config.h.in
  ${CMAKE_BINARY_DIR}/include/OpenMS/test_config.h
)

# Source files
set(DATE_SOURCES
  ../src/openms/source/DATASTRUCTURES/Date.cpp
)

# Create Date library
add_library(OpenMS_Date STATIC ${DATE_SOURCES})
target_link_libraries(OpenMS_Date PUBLIC Qt5::Core)
target_include_directories(OpenMS_Date PUBLIC ${CMAKE_BINARY_DIR}/include)

# Test executable
add_executable(Date_test ../src/tests/class_tests/openms/source/Date_test.cpp)
target_link_libraries(Date_test OpenMS_Date Qt5::Core)
target_include_directories(Date_test PUBLIC ${CMAKE_BINARY_DIR}/include)

# Enable testing
enable_testing()
add_test(NAME Date_test COMMAND Date_test)
