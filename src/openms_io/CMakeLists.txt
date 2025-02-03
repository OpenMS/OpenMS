project("OpenMS_IO")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#------------------------------------------------------------------------------
# includes for sources
include(${PROJECT_SOURCE_DIR}/includes.cmake)

#------------------------------------------------------------------------------
# Dependencies
set(OPENMS_IO_DEP_LIBRARIES
    OpenMS
    XercesC::XercesC
    Qt6::Core
    Qt6::Network
)

set(OPENMS_IO_DEP_PRIVATE_LIBRARIES
    BZip2::BZip2
    Boost::boost
    Boost::iostreams
    ZLIB::ZLIB
)

#------------------------------------------------------------------------------
# Create library
openms_add_library(TARGET_NAME  OpenMS_IO
                   SOURCE_FILES  ${OpenMS_IO_sources}
                   HEADER_FILES  ${OpenMS_IO_sources_h}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include
                                   ${CMAKE_CURRENT_BINARY_DIR}/include
                   LINK_LIBRARIES ${OPENMS_IO_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${OPENMS_IO_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "OpenMS_IO/")

if (MSVC)
  target_compile_options(OpenMS_IO PRIVATE "/we4100")
  target_compile_options(OpenMS_IO PRIVATE "/we4189")
endif()