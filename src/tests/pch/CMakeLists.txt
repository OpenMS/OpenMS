# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Timo Sachsenberg $

cmake_minimum_required(VERSION 3.24)
project(OpenMSPCHTest LANGUAGES C CXX)
enable_testing()

include(../../../cmake/precompiled_headers.cmake)

add_library(pch_probe STATIC probe.cpp probe.c)
target_compile_features(pch_probe PUBLIC cxx_std_23)
openms_add_precompiled_headers(pch_probe)

add_library(pch_disabled STATIC probe.cpp)
target_compile_features(pch_disabled PRIVATE cxx_std_23)
set_target_properties(pch_disabled PROPERTIES DISABLE_PRECOMPILE_HEADERS ON)
openms_add_precompiled_headers(pch_disabled)
get_target_property(_disabled_pch pch_disabled PRECOMPILE_HEADERS)
if(_disabled_pch)
  message(FATAL_ERROR "Target-level PCH opt-out was ignored")
endif()

# A dependent target must compile without inheriting the producer's PCH.
# The producer also exercises mixed C and C++ sources.
add_executable(pch_consumer main.cpp)
target_link_libraries(pch_consumer PRIVATE pch_probe)
add_test(NAME pch_consumer COMMAND pch_consumer)

get_target_property(_public_pch pch_probe INTERFACE_PRECOMPILE_HEADERS)
if(_public_pch)
  message(FATAL_ERROR "PCH headers leaked into the public interface")
endif()
get_target_property(_private_pch pch_probe PRECOMPILE_HEADERS)
if(OPENMS_USE_PCH AND NOT _private_pch)
  message(FATAL_ERROR "PCH was requested but not enabled")
elseif(NOT OPENMS_USE_PCH AND _private_pch)
  message(FATAL_ERROR "PCH was enabled by default")
endif()

install(TARGETS pch_probe EXPORT PCHProbeTargets ARCHIVE DESTINATION lib)
install(EXPORT PCHProbeTargets DESTINATION lib/cmake/PCHProbe)
