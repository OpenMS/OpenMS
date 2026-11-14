cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

### example CMakeLists.txt to develop programs using OpenMS
project("OpenMS_ExternalCodeTest")

include(CTest)

message("Received CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

## list all your executables here (a corresponding .cpp file should exist, e.g. TestExternalCode.cpp)
set(my_executables
	TestExternalCode
)

## list all classes here, which are required by your executables
## (all these classes will be linked into a library)
set(my_sources
	ExampleLibraryFile.cpp
)

## Find the OpenMS package. It provides the imported targets OpenMS::OpenMS and
## OpenMS::OpenSwathAlgo, OpenMS::OpenMS_CLI when the installation includes the TOPP
## tool framework (OpenMS_CLI_FOUND; a TOPP-style tool requests COMPONENTS CLI) and,
## with COMPONENTS GUI, OpenMS::OpenMS_GUI of the GUI layer (OpenMS-TOPP-GUI) when that
## is installed into the same prefix; its OpenMSGUI package also finds the Qt6 modules.
## Select an installation with -DOpenMS_DIR=<prefix>/lib/cmake/OpenMS
## (<prefix>/CMake on Windows) or by adding <prefix> to CMAKE_PREFIX_PATH.
find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI OPTIONAL_COMPONENTS GUI)
message(STATUS "Found OpenMS ${OpenMS_VERSION} at ${OpenMS_DIR}")

#------------------------------------------------------------------------------
## Package contract. This project doubles as the installed-consumer test of the
## OpenMS build (OPENMS_TEST_INSTALLED_CONSUMER) for the full development
## installation (core_only/ covers a core-only one), so it checks what the package
## promises before using it. A plain consumer project does not need this block.

## the CLI layer of the full installation was found ...
if(NOT OpenMS_CLI_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake found the required CLI component without setting OpenMS_CLI_FOUND")
endif()
## ... the libraries are namespaced imported targets ...
foreach(_target IN ITEMS OpenMS::OpenMS OpenMS::OpenSwathAlgo OpenMS::OpenMS_CLI)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "OpenMSConfig.cmake did not provide the target ${_target}")
  endif()
  get_target_property(_imported ${_target} IMPORTED)
  if(NOT _imported)
    message(FATAL_ERROR "${_target} is not an imported target")
  endif()
endforeach()
## ... the un-namespaced names of earlier releases remain as aliases ...
foreach(_alias IN ITEMS OpenMS OpenSwathAlgo OpenMS_CLI)
  if(NOT TARGET ${_alias})
    message(FATAL_ERROR "OpenMSConfig.cmake did not provide the compatibility alias ${_alias}")
  endif()
endforeach()
## ... and the bundled third-party libraries exist under the OpenMS:: namespace only
foreach(_leaked IN ITEMS IsoSpec Evergreen eol-bspline opentims_cpp openms_thermo_bridge)
  if(TARGET ${_leaked})
    message(FATAL_ERROR "OpenMSConfig.cmake leaked the un-namespaced third-party target ${_leaked}")
  endif()
endforeach()

## Qt6 is not part of the core library's link interface: a consumer that wants Qt finds
## and links it itself
get_target_property(_core_interface OpenMS::OpenMS INTERFACE_LINK_LIBRARIES)
if(_core_interface MATCHES "Qt6::")
  message(FATAL_ERROR "OpenMS::OpenMS leaks Qt6 into its link interface: ${_core_interface}")
endif()
## ... and neither is Eigen, for a shared installation: no installed header includes it,
## and CMake omits the PRIVATE dependencies of a SHARED library from the export, so a
## consumer never has to resolve Eigen3::Eigen. A STATIC installation legitimately keeps
## it as $<LINK_ONLY:Eigen3::Eigen>, so ask the imported target how it was built rather
## than trusting a recorded flag -- TYPE is what CMake itself reports, and it is also how
## a consuming project should answer this question. TestExternalCode_noeigen_configure
## configures this project with Eigen made unfindable, which is the other half of the check.
get_target_property(_core_type OpenMS::OpenMS TYPE)
if(_core_type STREQUAL "SHARED_LIBRARY" AND _core_interface MATCHES "Eigen")
  message(FATAL_ERROR "OpenMS::OpenMS leaks Eigen into its link interface: ${_core_interface}")
endif()
## ... and not through the include directories either. The link interface is not the only
## way Eigen can reach a consumer: a soft lookup reinstated in OpenMSConfig.cmake could
## append Eigen's include directories to the imported target, which leaves the link
## interface clean and still hands everyone a path <Eigen/Core> resolves from. Combined
## with an installed header including Eigen that compiles here and breaks for anyone
## without Eigen installed. Asked functionally rather than by matching path names, so a
## build directory that happens to contain "eigen" does not trip it.
## The package reports its header locations through a CMake file set, whose base
## directories reach INTERFACE_INCLUDE_DIRECTORIES as $<BUILD_INTERFACE:...>, so
## the raw property holds no directory to probe and this check would pass on
## every installation; openms_package_include_dirs() resolves them to paths.
include("${CMAKE_CURRENT_LIST_DIR}/openms_package_include_dirs.cmake")
openms_package_include_dirs(OpenMS::OpenMS _core_includes)
if(_core_type STREQUAL "SHARED_LIBRARY" AND _core_includes)
  foreach(_inc IN LISTS _core_includes)
    if(EXISTS "${_inc}/Eigen/Core")
      message(FATAL_ERROR "OpenMS::OpenMS exposes Eigen through its include directories: ${_inc}")
    endif()
  endforeach()
endif()

## the directory variables describe this installation, not the prefix of a dependency
if(NOT EXISTS "${OPENMS_DATA_DIR}/CHEMISTRY/unimod.xml")
  message(FATAL_ERROR "OPENMS_DATA_DIR='${OPENMS_DATA_DIR}' does not hold the OpenMS shared data")
endif()
if(NOT IS_DIRECTORY "${OPENMS_LIB_DIR}")
  message(FATAL_ERROR "OPENMS_LIB_DIR='${OPENMS_LIB_DIR}' is not a directory")
endif()
file(REAL_PATH "${OPENMS_LIB_DIR}" _lib_dir)
get_target_property(_configurations OpenMS::OpenMS IMPORTED_CONFIGURATIONS)
if(NOT _configurations)
  message(FATAL_ERROR "OpenMS::OpenMS has no imported configurations")
endif()
foreach(_configuration IN LISTS _configurations)
  get_target_property(_location OpenMS::OpenMS IMPORTED_LOCATION_${_configuration})
  get_filename_component(_location_dir "${_location}" DIRECTORY)
  file(REAL_PATH "${_location_dir}" _location_dir)
  if(NOT _location_dir STREQUAL _lib_dir)
    message(FATAL_ERROR "OPENMS_LIB_DIR='${OPENMS_LIB_DIR}' but the ${_configuration} OpenMS library is ${_location}")
  endif()
endforeach()

if(UNIX)
  set(_versioned_targets OpenMS::OpenMS OpenMS::OpenSwathAlgo OpenMS::OpenMS_CLI)
  foreach(_target IN LISTS _versioned_targets)
    get_target_property(_target_configs ${_target} IMPORTED_CONFIGURATIONS)
    foreach(_config IN LISTS _target_configs)
      get_target_property(_soname ${_target} IMPORTED_SONAME_${_config})
      string(REPLACE "OpenMS::" "" _name "${_target}")
      if(APPLE)
        set(_expected "@rpath/lib${_name}.${OpenMS_VERSION}.dylib")
      else()
        set(_expected "lib${_name}.so.${OpenMS_VERSION}")
      endif()
      if(NOT _soname STREQUAL _expected)
        message(FATAL_ERROR "${_target} SONAME '${_soname}' differs from '${_expected}'")
      endif()
    endforeach()
  endforeach()
endif()

## GUI component: OpenMSConfig.cmake delegates it to the OpenMSGUI package of the GUI
## layer, which OpenMS-TOPP-GUI installs next to it. Without that package (as in the
## installation of the installed-consumer test) it must not be found; the GUI layer has
## its own consumer test in OpenMS-TOPP-GUI.
if(NOT EXISTS "${OpenMS_DIR}/../OpenMSGUI/OpenMSGUIConfig.cmake"
   AND NOT EXISTS "${OpenMS_DIR}/OpenMSGUIConfig.cmake")
  if(OpenMS_GUI_FOUND OR OpenMS_WITH_GUI OR TARGET OpenMS::OpenMS_GUI)
    message(FATAL_ERROR "OpenMSConfig.cmake reports the GUI component for an installation without the GUI layer")
  endif()
elseif(OpenMS_GUI_FOUND AND NOT TARGET OpenMS::OpenMS_GUI)
  message(FATAL_ERROR "OpenMSConfig.cmake reports the GUI component without providing OpenMS::OpenMS_GUI")
endif()

#------------------------------------------------------------------------------
## library with additional classes from above
add_library(my_custom_lib STATIC ${my_sources})
## Projects written for earlier releases link the un-namespaced name. This library
## keeps that spelling on purpose, so the compatibility alias stays covered.
target_link_libraries(my_custom_lib OpenMS)
## disable auto-linking of boost (mainly needed for Windows. CMake takes care of choosing the right library)
target_compile_definitions(my_custom_lib PUBLIC BOOST_ALL_NO_LIB)

## add targets for the executables
foreach(i ${my_executables})
  # create the executable
  add_executable(${i} ${i}.cpp)
  ## link executables against OpenMS
  target_link_libraries(${i} OpenMS::OpenMS my_custom_lib)
endforeach(i)

## A TOPP-style tool of an external project: derives from TOPPBase, so it links the
## tool framework OpenMS::OpenMS_CLI (which carries OpenMS::OpenMS transitively).
add_executable(TestExternalCodeCLI TestExternalCodeCLI.cpp)
target_link_libraries(TestExternalCodeCLI OpenMS::OpenMS_CLI)
add_test(TestExternalCodeCLI TestExternalCodeCLI -greeting "hello from OpenMS_CLI")
set_tests_properties(TestExternalCodeCLI PROPERTIES PASS_REGULAR_EXPRESSION "hello from OpenMS_CLI")

## OSWPQ storage and the real OpenDIA source must build against installed public
## headers, with no Arrow target or source-tree include directory.
add_executable(TestExternalCodeOSWParquet TestExternalCodeOSWParquet.cpp)
target_link_libraries(TestExternalCodeOSWParquet OpenMS::OpenMS)
add_test(TestExternalCodeOSWParquet TestExternalCodeOSWParquet)
add_executable(TestExternalOpenDIA ../../topp/OpenDIA.cpp)
target_link_libraries(TestExternalOpenDIA OpenMS::OpenMS_CLI)
add_test(TestExternalOpenDIA TestExternalOpenDIA -test -write_ini "${CMAKE_CURRENT_BINARY_DIR}/external_opendia.ini")

## Enable testing - for Nightly Build log
add_test(TestExternalCode TestExternalCode)

if(UNIX AND BUILD_TESTING)
  # This unversioned, incompatible library simulates another OpenMS installation
  # found first via LD_LIBRARY_PATH / DYLD_LIBRARY_PATH.
  add_library(ConflictingOpenMS SHARED ConflictingOpenMS.cpp)
  set_target_properties(ConflictingOpenMS PROPERTIES
    OUTPUT_NAME OpenMS
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/conflicting-libraries")
  add_dependencies(TestExternalCode ConflictingOpenMS)
  if(APPLE)
    set(_library_path_variable DYLD_LIBRARY_PATH)
  else()
    set(_library_path_variable LD_LIBRARY_PATH)
    target_link_options(TestExternalCode PRIVATE "LINKER:--enable-new-dtags")
  endif()
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/soname-fallback")
  add_test(NAME TestExternalCode_versioned_soname_fallback
    COMMAND "$<TARGET_FILE:TestExternalCode>")
  set_tests_properties(TestExternalCode_versioned_soname_fallback PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/soname-fallback"
    ENVIRONMENT_MODIFICATION "${_library_path_variable}=path_list_prepend:$<TARGET_FILE_DIR:ConflictingOpenMS>")
endif()

## A WITH_ONNX installation ships the PeptDeep inference classes and their ONNX
## models under OPENMS_DATA_DIR/models, and the ONNX Runtime library among its
## dependencies. The consumer loads every model from the installation and runs a
## prediction, so a package that carries the classes but not the models, or a
## runtime that cannot be loaded, fails here rather than at a user's.
if(OpenMS_WITH_ONNX)
  set(_peptdeep_models_dir "${OPENMS_DATA_DIR}/models")
  foreach(_model IN ITEMS peptdeep_rt_dynamic.onnx peptdeep_ccs_dynamic.onnx peptdeep_ms2_dynamic.onnx)
    if(NOT EXISTS "${_peptdeep_models_dir}/${_model}")
      message(FATAL_ERROR "OpenMS_WITH_ONNX is set but the installation does not ship ${_peptdeep_models_dir}/${_model}")
    endif()
  endforeach()
  add_executable(TestExternalCodePeptDeep TestExternalCodePeptDeep.cpp)
  target_link_libraries(TestExternalCodePeptDeep OpenMS::OpenMS)
  add_test(TestExternalCodePeptDeep TestExternalCodePeptDeep "${_peptdeep_models_dir}")
endif()

## Windows has no RPATH: executables find OpenMS.dll (and the DLLs of the
## dependencies OpenMS was built against) through PATH. OpenMSConfig.cmake
## tells us where the OpenMS libraries were installed.
if(WIN32)
  foreach(_test IN ITEMS TestExternalCode TestExternalCodeCLI TestExternalCodePeptDeep TestExternalCodeOSWParquet TestExternalOpenDIA)
    if(TEST ${_test})
      set_property(TEST ${_test} APPEND PROPERTY
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${OPENMS_LIB_DIR}")
    endif()
  endforeach()
endif()
