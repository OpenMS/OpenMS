# CMake installation file
# Created 2008/10 during OpenMS retreat by CB && AB

PROJECT("OpenMS")

########################################################
### entries meant to be configured using CMake cache ###
### - do NOT hardcode them here!                     ###
### - edit them within CMakeCache.txt using ccmake	 ###
########################################################
# CMAKE_FIND_ROOT_PATH
# CMAKE_BUILD_TYPE
# STL_DEBUG
# DB_TEST
# QT_DB_PLUGIN
# MT_CUDA_BUILD_TYPE

########################################################
###    manual entries (edit this for new release)    ###
########################################################

set(OPENMS_PACKAGE_VERSION_MAJOR "1")
set(OPENMS_PACKAGE_VERSION_MINOR "9")
set(OPENMS_PACKAGE_VERSION_PATCH "0")

########################################################
###            end manual entries                    ###
########################################################

## Heart of the BUILD system : only edit when you know what you're doing (we don't)
## quick manual for most commands: http://www.cmake.org/cmake/help
## useful predefined variables: http://www.paraview.org/Wiki/CMake_Useful_Variables
##
## naming conventions:
##
## prefix a variable with 'CF_' if it is used to configure a file!
## e.g., CF_LibOpenMSExport


cmake_minimum_required(VERSION 2.6.3 FATAL_ERROR)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set(CF_OPENMS_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}" CACHE INTERNAL "OpenMS VERSION" FORCE)

## SOME MAC OS X STUFF
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # force cmake > 2.8.1 to avoid bugs with OSX_DEPLOYMENT_TARGET
  cmake_minimum_required(VERSION 2.8.1 FATAL_ERROR)
	# due to some awkward updates of cmake in the 2.8.x series we need to add the following
  find_library(CARBON NAMES Carbon)
  find_package(OpenGL)
  set(APPLE_EXTRA_LIBS ${CARBON} ${OPENGL_LIBRARY})
endif()


## sanity check: sometimes CMAKE_SIZEOF_VOID_P just vanishes (when updating CMake).
if (NOT CMAKE_SIZEOF_VOID_P)
	message(FATAL_ERROR "'CMAKE_SIZEOF_VOID_P' is undefined. Thus you should delete CMakeFiles (the directory) and the CMakeCache.txt and rerun CMake again! This is some weird CMake bug that seems to appear when updating the CMake version.")
endif()

if (CMAKE_SIZEOF_VOID_P MATCHES "8")
	set(OPENMS_64BIT_ARCHITECTURE 1 CACHE INTERNAL "Architecture-bits")
	message(STATUS "Architecture: 64 bit")
else()
	set(OPENMS_64BIT_ARCHITECTURE 0 CACHE INTERNAL "Architecture-bits")
	message(STATUS "Architecture: 32 bit")
endif()

# some of our own macros (OPENMS_CHECKLIB, QT4_WRAP_UI_OWN, ...)
include (cmake/OpenMSBuildSystem_macros.cmake)


## Set default build type (if not set by user on command line)
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
## force build type into the cache (needs to be set beforehand)
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
	set(CF_OPENMS_ASSERTIONS 1)
else()
	set(CF_OPENMS_ASSERTIONS 0)
endif()
set(CF_OPENMS_ASSERTIONS ${CF_OPENMS_ASSERTIONS} CACHE INTERNAL "Enables debug messages (precondition and postconditions are enabled, a bit slower) - this is NOT changing any compiler flags!" FORCE)


########################################################
###    compiler flags																 ###
########################################################

## Fill this with compile flags that external projects should use as well
## for OpenMS internal flags (not promoted to external compiler flags) append to CMAKE_CXX_FLAGS
set(CF_OPENMS_ADDCXX_FLAGS) ## see OpenMS/cmake/OpenMSConfig.cmake.in to see how its configured and used (i.e. as OPENMS_ADDCXX_FLAGS)

if (CMAKE_COMPILER_IS_GNUCXX)

	# determine gcc version
	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GNUCXX_VERSION)
	string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GNUCXX_VERSION})
	list(GET GCC_VERSION_COMPONENTS 0 GNUCXX_MAJOR_VERSION)
	list(GET GCC_VERSION_COMPONENTS 1 GNUCXX_MINOR_VERSION)

  add_definitions(-Wall -Wextra -Wno-non-virtual-dtor -Wno-long-long -Wno-variadic-macros)
  if (NOT MT_ENABLE_CUDA)  # necessary since CUDA contains non-pedantic code
		add_definitions(--pedantic)
	endif()

	## Recommended setting for eclipse, see http://www.cmake.org/Wiki/CMake:Eclipse
	if (CMAKE_GENERATOR STREQUAL "Eclipse CDT4 - Unix Makefiles")
		add_definitions(-fmessage-length=0)
	endif()

	if (NOT OPENMS_64BIT_ARCHITECTURE AND ${GNUCXX_MAJOR_VERSION} MATCHES "4" AND ${GNUCXX_MINOR_VERSION} MATCHES "3")
		add_definitions(-march=i486)
	endif()

elseif (MSVC)
	## do not use add_definitions
	## add definitions also lands in stuff like RC_DEFINITION which tend to fail if you use
	## Eclipse CDT 4 - NMAKE generator
	## use set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} ...") instead

	## disable dll-interface warning
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4251 /wd4275")

	## disable deprecated functions warning (e.g. for POSIX functions)
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4996")

	## disable explicit template instantiation request for partially defined classes
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4661")

	## disable warning: "decorated name length exceeded, name was truncated"
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4503")
	
	## don't warn about unchecked std::copy()
	add_definitions(/D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_DEPRECATE)

	## xerces bug workaround
	add_definitions(/DOPENMS_XERCESDLL)

	## FeatureFinder.obj is huge and won't compile in VS2008 debug otherwise:
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")

	## use multiple Cores (if available)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

  if (NOT OPENMS_64BIT_ARCHITECTURE)
    ## enable SSE1 on 32bit, on 64bit the compiler flag does not exist
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE")
  endif()
elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_COMPILER_IS_CLANG true CACHE INTERNAL "Is CLang compiler (clang++)")
else()
	set(CMAKE_COMPILER_IS_INTELCXX true CACHE INTERNAL "Is Intel C++ compiler (icpc)")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CF_OPENMS_ADDCXX_FLAGS}")

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
	CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
	if (WITH_FPIC)
		add_definitions(-fPIC)
	endif()
endif()

## -Wconversion flag for GCC
set(CXX_WARN_CONVERSION OFF CACHE BOOL "Enables warnings for type conversion problems (GCC only)")
if (CXX_WARN_CONVERSION)
	if (CMAKE_COMPILER_IS_GNUCXX)
		add_definitions(-Wconversion)
	endif()
endif()
Message(STATUS "Compiler checks for conversion: ${CXX_WARN_CONVERSION}")

## STL-DEBUG (only for GCC and in debug mode)
set(STL_DEBUG OFF CACHE BOOL "[GCC only] Enable STL-DEBUG mode (very slow).")
if (STL_DEBUG)
  if (CMAKE_COMPILER_IS_GNUCXX)
		if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
			# add compiler flag
    	add_definitions(/D_GLIBCXX_DEBUG)
    	Message(STATUS "STL debug mode: ${STL_DEBUG}")
	  else()
	    Message(WARNING "STL debug mode is supported for OpenMS debug mode only")
	  endif()
  else()
    Message(WARNING "STL debug mode is supported for compiler GCC only")
  endif()
else()
	Message(STATUS "[GCC only] STL debug mode: ${STL_DEBUG}")
endif()


#########################################################
#							 Useful third party libs (optional)       #
#########################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

## CUDA checking
set(MT_ENABLE_CUDA OFF CACHE BOOL "Should CUDA support be enabled (version <= 2.1 currently supported).")
set(MT_CUDA_BUILD_TYPE Device CACHE STRING "Switch between Device and Emulation mode.")
if (MT_ENABLE_CUDA)
	MESSAGE(STATUS "NVIDIA CUDA: ${MT_CUDA_BUILD_TYPE} mode")
	find_package(Cuda)
else()
  MESSAGE(STATUS "NVIDIA CUDA: ${MT_ENABLE_CUDA}")
endif()

if (FOUND_CUDART AND MT_ENABLE_CUDA) #CUDART and not CUDA!!!
	add_definitions(/DOPENMS_HAS_CUDA)
	MESSAGE (STATUS "Found CUDA header files in: ${FOUND_CUDA_NVCC_INCLUDE}") #Not CUDA_INCLUDE_DIR
	MESSAGE (STATUS "Found CUDA library at: ${FOUND_CUDART}")
	INCLUDE_DIRECTORIES(${FOUND_CUDA_NVCC_INCLUDE})
  set (CUDA_NVCC_INCLUDE_ARGS ${CUDA_NVCC_INCLUDE_ARGS} -I ${PROJECT_BINARY_DIR}/include/)
endif()

## TBB
set(MT_ENABLE_TBB OFF CACHE BOOL "Should Intel Threading Building Blocks support be enabled.")
set(MT_TBB_INCLUDE_DIR CACHE PATH "Intel Threading Building Blocks 'include' directory.")
set(MT_TBB_LIBRARY_DIR CACHE PATH "Intel Threading Building Blocks libraries directory.")
MESSAGE(STATUS "Intel TBB: ${MT_ENABLE_TBB}")
if (MT_ENABLE_TBB)
	find_package(TBB)
endif()

if (TBB_FOUND)
	INCLUDE_DIRECTORIES(${MT_TBB_INCLUDE_DIR})
	add_definitions(/DOPENMS_HAS_TBB)
endif()

## OpenMP
set(MT_ENABLE_OPENMP ON CACHE BOOL "Should OpenMP support be enabled")
if (MT_ENABLE_OPENMP)
	find_package(OpenMP)
endif()
MESSAGE(STATUS "OPENMP: ${MT_ENABLE_OPENMP}")

if (OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}") ## do NOT use add_definitions() here, because RC.exe on windows will fail
	if (NOT MSVC)
		set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()

########################################################
###    external libs (contrib or system)						 ###
########################################################

include(cmake/OpenMSBuildSystem_externalLibs.cmake)


########################################################
###    useful programms															 ###
########################################################

find_program(SVNVERSION_EXECUTABLE
             svnversion
             PATHS "c:/programme/subversion/bin" "c:/program files/subversion/bin"   ### additional search paths (along with $PATH)
             DOC "svnversion executable which helps in determining the svn revision when building TOPP tools")

if ("${SVNVERSION_EXECUTABLE}" STREQUAL "SVNVERSION_EXECUTABLE-NOTFOUND")
	message(STATUS "Info: The programm svnversion could not be found. SVN-revision information will not be available! Add the location of svnversion(.exe) to your PATH environment variable if you require SVN-revision.")
	set(OPENMS_HAS_SVNVERSION OFF CACHE INTERNAL "SVNVersion(.exe) present?")
else()
	set(OPENMS_HAS_SVNVERSION ON CACHE INTERNAL "SVNVersion(.exe) present?")
endif()

set(SVN_REVISION_FILE "${PROJECT_BINARY_DIR}/include/OpenMS/openms_svn_revision.h")



########################################################
###    configure config.h														 ###
########################################################

## this need to be set before config.h is configured!
set(CMAKE_INSTALL_PREFIX "" CACHE INTERNAL "This must be equal to INSTALL_PREFIX!" FORCE)
set(INSTALL_PREFIX "." CACHE PATH "Installation prefix for OpenMS and TOPP")
## are we building a shared or static lib?! (BOTH within the same BUILD-tree is NOT possible with OpenMS!!)
set(BUILD_SHARED_LIBS true)

include(cmake/OpenMSBuildSystem_configh.cmake)

## add configured files (used in 'cmake/includes.cmake')
set (OPENMS_CONFIGURED_FILES ${SVN_REVISION_FILE} ${CONFIGURED_CONFIG_H} ${CONFIGURED_OPENMS_PACKAGE_VERSION_H})

## Configures what is tested (DB-test, MS2-search engines)
include(cmake/OpenMSBuildSystem_testConfig.cmake)



########################################################
###    Documentation                                 ###
########################################################

include(cmake/OpenMSBuildSystem_doc.cmake)



########################################################
###    BUILD the lib                                 ###
########################################################

#package type
set(VALID_PACKAGE_TYPES "none" "rpm" "deb" "dmg")
set(PACKAGE_TYPE "none" CACHE STRING "Package type (internal): ${VALID_PACKAGE_TYPES}")
# check build type
list(FIND VALID_PACKAGE_TYPES ${PACKAGE_TYPE} list_pos)
if( ${list_pos} EQUAL -1 )
	message(STATUS "The PACKAGE_TYPE ${PACKAGE_TYPE} is invalid")
	message(STATUS "Valid PACKAGE_TYPEs are:")
	foreach( BT ${VALID_PACKAGE_TYPES} )
		message(STATUS " * ${BT}")
	endforeach()
	message(FATAL_ERROR "Aborting ...")
endif()


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(OPENMS_BINARY_DIR "${PROJECT_BINARY_DIR}/bin")
set(OPENMS_WIN32_DLL_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

## big include file for headers and C files, which fills the OpenMS_sources variable
include (cmake/includes.cmake)

if (MSVC)
	## use OpenMSd.dll in debug mode
	SET(CMAKE_DEBUG_POSTFIX d)
endif()

## add library target
## warning: set BUILD_SHARED_LIBS to decide if library is shared or static (see above)! We need the BUILD_SHARED_LIBS flag to set declspec flags for MSVC!
if (FOUND_CUDART AND MT_ENABLE_CUDA) # we need both conditions due to possible present cached entries
	CUDA_ADD_LIBRARY(OpenMS ${OpenMS_sources} ${Cuda_sources})
else()
	ADD_LIBRARY(OpenMS ${OpenMS_sources})
endif()

### OpenMS_GUI lib
ADD_LIBRARY(OpenMS_GUI ${OpenMSVisual_sources})
if (MSVC)
	set (GUI_lnk_flags "/FORCE:MULTIPLE /INCREMENTAL:NO /ignore:4006 /ignore:4088")
	set_target_properties(OpenMS_GUI PROPERTIES LINK_FLAGS ${GUI_lnk_flags}) ## allow multiple definitions of symbols (e.g. from template instanciations or STL-derived classes)
endif()	

## update the openms_svn_revision.h file before compiling OpenMS
add_custom_target(svn_revision_update ALL
	    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -DSVN_REVISION_FILE=${SVN_REVISION_FILE} -DOPENMS_HAS_SVNVERSION=${OPENMS_HAS_SVNVERSION} -DSVNVERSION_EXECUTABLE=${SVNVERSION_EXECUTABLE}
	                             -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/svnversion.cmake)
add_dependencies(OpenMS svn_revision_update)
set_source_files_properties(${SVN_REVISION_FILE} PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE)

## all the dependency libraries are linked into libOpenMS.so, except Qt and CUDA which are still dynamic
set(OPENMS_DEP_LIBRARIES ${CONTRIB_CBC} ${CONTRIB_GSL} ${CONTRIB_GSLCBLAS} ${DEP_LIBSVM_LIBRARY} ${CONTRIB_XERCESC} ${Boost_LIBRARIES} ${CONTRIB_BZIP2} ${CONTRIB_ZLIB} ${CONTRIB_GLPK})
if (TBB_FOUND)
	list(APPEND OPENMS_DEP_LIBRARIES ${TBB_LIBRARIES})
endif()
if (MSVC)
	list(APPEND OPENMS_DEP_LIBRARIES opengl32.lib)
endif()

target_link_libraries(OpenMS ${APPLE_EXTRA_LIBS} ${QT_LIBRARIES} ${OPENMS_DEP_LIBRARIES})
target_link_libraries(OpenMS_GUI OpenMS ${QT_LIBRARIES})

## OPENMS_LIBRARIES defines the libraries used by OpenMS; this should be used to link against executables
## somehow the link dependencies of cmake are broken such that when using POSTFIX names for libs the dependencies are not forwarded...
## we fix this by adding the dependencies of OpenMS directly to the executables as well.
set(OPENMS_LIBRARIES optimized OpenMS${CMAKE_RELEASE_POSTFIX} debug OpenMS${CMAKE_DEBUG_POSTFIX} ${QT_LIBRARIES})
if (MSVC) ## Windows
	list(APPEND OPENMS_LIBRARIES ${OPENMS_DEP_LIBRARIES})
endif()

## link libraries for GUI Tools
set(OPENMS_GUI_LIBRARIES optimized OpenMS_GUI${CMAKE_RELEASE_POSTFIX} debug OpenMS_GUI${CMAKE_DEBUG_POSTFIX} OpenMS ${QT_LIBRARIES})

## directory for OpenMS(d).lib (required for linking executables)
if (MSVC)
	link_directories(${PROJECT_BINARY_DIR})
endif()

####### TOPP #########
set(TOPP_DIR source/APPLICATIONS/TOPP)
include(${TOPP_DIR}/executables.cmake)
foreach(i ${TOPP_executables})
	add_executable(${i} ${TOPP_DIR}/${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
	if (OPENMP_FOUND AND NOT MSVC)
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()
endforeach(i)
add_custom_target(TOPP)
add_dependencies(TOPP ${TOPP_executables})

##### UTILS ########
set(UTILS_DIR source/APPLICATIONS/UTILS)
set(UTILS_executables)
include(${UTILS_DIR}/executables.cmake)
foreach(i ${UTILS_executables})
	add_executable(${i} ${UTILS_DIR}/${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(UTILS)
add_dependencies(UTILS ${UTILS_executables})

##### GUI tools ########
set(GUI_DIR source/VISUAL/APPLICATIONS/GUITOOLS)
include(${GUI_DIR}/executables.cmake)
foreach(i ${GUI_executables})
	set(resource_file ${CMAKE_SOURCE_DIR}/${GUI_DIR}/${i}.rc)
	set(resource_dir ${CMAKE_SOURCE_DIR}/${GUI_DIR}/${i}-resources/)
	## add icons to TOPPView and INIFileEditor
	if (MSVC AND EXISTS ${resource_file})
		Message(STATUS "Setting resource file ${resource_file} for ${i}")
		add_executable(${i} ${GUI_DIR}/${i}.C ${resource_file})
	elseif(APPLE AND EXISTS ${resource_dir})
		add_executable(${i} ${GUI_DIR}/${i}.C)
		# extra work for pseudo-mac-bundle
		# move prepared package to RUNTIME_OUTPUT
		file(COPY "${PROJECT_SOURCE_DIR}/cmake/MacOSX/${i}-resources/${i}.app"
		DESTINATION "${PROJECT_BINARY_DIR}/bin")
	else()
		add_executable(${i} ${GUI_DIR}/${i}.C)
	endif()

	## append visual lib as dependency for GUI tools
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	if (OPENMP_FOUND AND NOT MSVC)
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()
endforeach(i)
add_custom_target(GUI)
add_dependencies(GUI ${GUI_executables})

## some regular TOPP/UTILS need the GUI lib
set(executables_withGUIlibDep ${TOPP_executables_with_GUIlib} ${UTILS_executables_with_GUIlib})
foreach(i ${executables_withGUIlibDep}) 
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
endforeach(i)

	
##### Doc progs ####
set(OpenMS_doc_executables)
include(doc/doxygen/parameters/executables.cmake)
foreach(i ${OpenMS_doc_executables})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "doc/doxygen/parameters")
	add_executable(${i} EXCLUDE_FROM_ALL doc/doxygen/parameters/${i}.C)
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")
endforeach(i)
add_custom_target(doc_progs)
add_dependencies(doc_progs ${OpenMS_doc_executables} TOPP UTILS GUI)

##### Examples #####
if ("${PACKAGE_TYPE}" STREQUAL "none")
	set(EXAMPLES_executables)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
	include(source/EXAMPLES/executables.cmake)
	foreach(i ${EXAMPLES_executables})
	  add_executable(${i} source/EXAMPLES/${i}.C)
	  target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	endforeach(i)
	add_custom_target(Tutorials_build)
	add_dependencies(Tutorials_build ${EXAMPLES_executables})

	if (MSVC_IDE)
		set(EXAMPLES_EXE_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$(OutDir))
	else()
		set(EXAMPLES_EXE_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
	endif()


	## Create single target for each example, add it as a dependency to the main target
	foreach(i ${EXAMPLES_executables})
		list(APPEND EX_LIST COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/source/EXAMPLES ${EXAMPLES_EXE_DIR}/${i})
	endforeach(i)
	list(APPEND EX_LIST "echo all done")
	add_custom_target(Tutorials_exec
		COMMENT "Executes the examples"
	)
	separate_arguments(EX_LIST)
	add_custom_command(TARGET Tutorials_exec
                    POST_BUILD
										${EX_LIST}
										)

endif()

########################################################
###                      Tests                       ###
########################################################
if("${PACKAGE_TYPE}" STREQUAL "none")
	INCLUDE(Dart) ## for Nighlty Build log
	add_subdirectory(source/TEST EXCLUDE_FROM_ALL) ## configure the Test project (but do not add it to OpenMS itself)

  ## CppCheck testing
  find_package(cppcheck)
  if( CPPCHECK_FOUND )
    include( cmake/cppcheck.cmake )
  endif()
endif()

########################################################
###    verbose Post-build messages and help targets  ###
########################################################

include(cmake/OpenMSBuildSystem_messages.cmake)

######################################################
#### build settings and configs for external code ####
######################################################
if ("${PACKAGE_TYPE}" STREQUAL "none")
	## grep definitions (-D<XXX>) that were set using 'add_definitions()'
	get_property(OPENMS_CMP_DEFS DIRECTORY . PROPERTY COMPILE_DEFINITIONS)
	## convert to CXX flags
	foreach(i ${OPENMS_CMP_DEFS})
		set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} -D${i}")
	endforeach(i)
	## for add_definition commands (currently empty, but who knows...)
	set(CF_OPENMS_ADD_DEFINITIONS)
	## the filename for library settings (this file will be included in OpenMSConfig.cmake)
	set(CF_LibOpenMSExport ${PROJECT_BINARY_DIR}/cmake/OpenMSLibraryExport.cmake)

	## create code that allows to call find_package()
	## (must come after the above code that sets the variables!)
	get_directory_property(OPENMS_INCLUDE_DIRS  INCLUDE_DIRECTORIES)
	configure_file(
		"${PROJECT_SOURCE_DIR}/cmake/OpenMSConfig.cmake.in"
		"${PROJECT_BINARY_DIR}/cmake/OpenMSConfig.cmake"
		@ONLY
	)

	## Create OpenMSLibraryExport.cmake:
	## Note: export_build_settings and
	##       export_library_dependencies are deprecated and should not be used
	if ("${INSTALL_PREFIX}" STREQUAL ".")
		set(CVERSION "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}")
		if (CVERSION VERSION_LESS "2.8.0")
			message(STATUS "Info: Cannot register OpenMS with CMake! For external code, set the path to OpenMS during find_package() manually.")
		else()
			## register OpenMS with CMake so that it can be found easily
			export(PACKAGE OpenMS)
		endif()
		## export our lib
		export(TARGETS OpenMS FILE ${CF_LibOpenMSExport})
	else()
		include(cmake/OpenMSBuildSystem_install.cmake)
	endif()

else()
######################################################
####          install/copy (cpack mainly)         ####
######################################################
	## general definitions
	set(CPACK_PACKAGE_NAME "OpenMS")
	set(CPACK_PACKAGE_VENDOR "OpenMS.de")
	set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMS - A framework for mass spectrometry")
	set(CPACK_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}")
	set(CPACK_PACKAGE_VERSION_MAJOR "${OPENMS_PACKAGE_VERSION_MAJOR}")
	set(CPACK_PACKAGE_VERSION_MINOR "${OPENMS_PACKAGE_VERSION_MINOR}")
	set(CPACK_PACKAGE_VERSION_PATCH "${OPENMS_PACKAGE_VERSION_PATCH}")
	set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMS-${CPACK_PACKAGE_VERSION}")
	set(CPACK_PACKAGE_DESCRIPTION_FILE ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageDescriptionFile.cmake)
	set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/License.txt)
	set(CPACK_RESOURCE_FILE_WELCOME ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceWelcomeFile.txt)
	set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceReadme.txt)

	# install routines for MacOSX
	if("${PACKAGE_TYPE}" STREQUAL "dmg")
		include(cmake/package_dragndrop_dmg.cmake)
	endif()

	if("${PACKAGE_TYPE}" STREQUAL "rpm")
		include(cmake/package_rpm.cmake)
	endif()

	if("${PACKAGE_TYPE}" STREQUAL "deb")
		include(cmake/package_deb.cmake)
	endif()

endif()

