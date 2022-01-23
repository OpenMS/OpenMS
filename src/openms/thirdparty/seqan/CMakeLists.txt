project("OpenMS-seqan")
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# export include directory
add_library(SeqAn INTERFACE)
target_include_directories(SeqAn SYSTEM INTERFACE 
                            "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>"
                            "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/thirdparty/seqan/include>")

install(TARGETS SeqAn EXPORT SeqAnTargets PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}/thirdparty/seqan COMPONENT Development)
install_directory(${CMAKE_CURRENT_LIST_DIR} ${INSTALL_INCLUDE_DIR}/thirdparty/seqan include)
# Exporting is done in the higher-level cmake files
#install(EXPORT SeqAnExport DESTINATION lib/cmake/seqan)
