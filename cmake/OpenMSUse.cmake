# - Use Module for external projects using OpenMS
# Sets up C and C++ to use OpenMS.  It is assumed that find_package(OpenMS)
# has already been loaded.

include(CMakeImportBuildSettings)

cmake_import_build_settings(${OpenMS_BUILD_SETTINGS_FILE})


# this is a hack (????), because build settings (compiler flags and Co.
# seem not to be automatically transferred to projects using the 
# build settings?!?
set(CMAKE_CXX_FLAGS ${CMAKE_BUILD_SETTING_CXX_FLAGS})
set(CMAKE_C_FLAGS ${CMAKE_BUILD_SETTING_C_FLAGS})



include_directories(${OPENMS_INCLUDE_DIRS})
link_directories(${OPENMS_LIBRARIES_DIRS})
