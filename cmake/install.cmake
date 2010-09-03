	
set(INSTALL_FORCE_DOC FALSE) ## do not force documentation to be present
include(cmake/install_common.cmake)

## create script that allows external projects to use our OpenMS lib
get_filename_component(LibOpenMSExportName "${CF_LibOpenMSExport}" NAME)
install(EXPORT OpenMSLibExportGroup DESTINATION cmake/ FILE ${LibOpenMSExportName})

