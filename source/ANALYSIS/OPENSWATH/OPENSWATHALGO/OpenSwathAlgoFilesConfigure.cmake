# configure OpenSwathAlgoConfig.h.in file for dll_export macros etc.
## replace any variables in config.h.in with current values
set (CONFIGURED_OPENSWATHALGOCONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h.in ${CONFIGURED_OPENSWATHALGOCONFIG_H})
