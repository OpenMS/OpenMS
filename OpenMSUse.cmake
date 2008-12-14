### use file for external projects
include(CMakeImportBuildSettings)

cmake_import_build_settings(${OpenMS_BUILD_SETTINGS_FILE})

include_directories(${OPENMS_INCLUDE_DIRS})

link_directories(${OPENMS_LIBRARIES_DIRS})
