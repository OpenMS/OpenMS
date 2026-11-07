# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(OpenDIA)

if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
endif()

add_library(OpenDIA_impl STATIC
  source/OpenDIASearch.cpp
  source/OpenDIARerank.cpp
  source/OpenDIATraces.cpp
  include/OpenDIASearch.h
  include/OpenDIARerank.h
  include/OpenDIATraces.h)
target_include_directories(OpenDIA_impl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(OpenDIA_impl PUBLIC OpenMS::OpenMS_CLI)
if(COMMAND openms_add_library_compiler_flags)
  openms_add_library_compiler_flags(OpenDIA_impl)
endif()

add_executable(OpenDIA OpenDIA.cpp)
target_link_libraries(OpenDIA PRIVATE OpenDIA_impl)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(OpenDIA)
endif()

add_executable(OpenDIARerank OpenDIARerank.cpp)
target_link_libraries(OpenDIARerank PRIVATE OpenDIA_impl)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(OpenDIARerank)
endif()

add_executable(OpenDIATraces OpenDIATraces.cpp)
target_link_libraries(OpenDIATraces PRIVATE OpenDIA_impl)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(OpenDIATraces)
endif()

option(OPENDIA_WITH_TORCH "Build experimental native OpenDIACNN training and scoring (requires CPU LibTorch)" OFF)
if(OPENDIA_WITH_TORCH)
  find_package(Torch CONFIG REQUIRED)
  # The standard-only backend API shares std::string with OpenMS, so its ABI must match.
  if(TORCH_CXX_FLAGS MATCHES "_GLIBCXX_USE_CXX11_ABI=0")
    message(FATAL_ERROR "OpenDIACNN requires a C++11-ABI LibTorch distribution compatible with OpenMS")
  endif()
  separate_arguments(_opendia_torch_flags NATIVE_COMMAND "${TORCH_CXX_FLAGS}")
  add_library(OpenDIA_cnn STATIC source/OpenDIACNN.cpp include/OpenDIACNN.h
    source/OpenDIACNNQValues.cpp include/OpenDIACNNQValues.h)
  target_compile_features(OpenDIA_cnn PUBLIC cxx_std_17)
  target_include_directories(OpenDIA_cnn PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_compile_options(OpenDIA_cnn PRIVATE ${_opendia_torch_flags})
  target_link_libraries(OpenDIA_cnn PUBLIC ${TORCH_LIBRARIES})
  add_executable(OpenDIACNN OpenDIACNN.cpp)
  target_link_libraries(OpenDIACNN PRIVATE OpenDIA_impl OpenDIA_cnn)
  if(COMMAND openms_add_executable_compiler_flags)
    openms_add_executable_compiler_flags(OpenDIACNN)
  endif()
endif()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(CTest)
  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()
endif()
