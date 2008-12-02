dnl -*- Mode: M4; tab-width: 2; -*-
dnl vi: set ts=2:
dnl
dnl    $Id: configure.ac,v 1.21 2006/05/30 15:46:38 marc_sturm Exp $
dnl
dnl    Process this file with autoconf to produce a configure script.
dnl   

dnl 
dnl    Some initialization stuff...
dnl 
AC_INIT(OpenMS, 1.2, www.OpenMS.de)
AC_CONFIG_HEADER([
	config.h:config/config.h.in
])
AC_LANG(C++)
PROJECT=OPENMS
define(PROJECTUPPER, OPENMS)
define(PROJECT, OpenMS)

dnl		
dnl    Define command line arguments for compiler,
dnl    compiler flags, includes and libraries.
dnl
AC_ARG_ENABLE(optimization,    [  --disable-optimization  disables optimizing the library for speed.])
AC_ARG_ENABLE(debuginfo,       [  --enable-debuginfo      add -g to the compiler flags (to create debug information). Disables optimization!])
AC_ARG_ENABLE(stl-debug,       [  --enable-stl-debug      enables debug mode of the STL when using gcc.])
AC_ARG_ENABLE(pedantic,        [  --disable-pedantic      removes the '-pedantic' flag when using gcc.])
AC_ARG_ENABLE(externalmemory,  [  --enable-externalmemory enable external memory structure for all MSExperiments (warning: performance hit). Only recommended if you run out of RAM+swap!])
AC_ARG_ENABLE(openmp,          [  --enable-openmp         enable usage of OpenMP features, which enables compiling OpenMS to use more than one thread.])
AC_ARG_ENABLE(db-plugin,       [  --with-db-plugin        set the Qt SQL DB plugin to use (default: QMYSQL)])
AC_ARG_WITH(compiler,          [  --with-compiler=CXX     use CXX as C++ compiler])
AC_ARG_WITH(cxxflags,          [  --with-cxxflags=FLAGS   add FLAGS to the C++ compiler flags (commas are converted to blanks)])
AC_ARG_WITH(ldflags,           [  --with-ldflags=FLAGS    add FLAGS to the linker flags (commas are converted to blanks)])
AC_ARG_WITH(arflags,           [  --with-arflags=FLAGS    add FLAGS to the flags for the creation of the static libraries])
AC_ARG_WITH(dynarflags,        [  --with-dynarflags=FLAGS add FLAGS to the flags for the creation of the shared libraries])
AC_ARG_WITH(qt,                [  --with-qt=QTDIR         Qt4 is installed in QTDIR. Prefer this over the following four options.])
AC_ARG_WITH(qt-incl,           [  --with-qt-incl=DIR      Qt4 header files are in DIR])
AC_ARG_WITH(qt-libs,           [  --with-qt-libs=DIR      Qt4 libraries are in DIR])
AC_ARG_WITH(moc,               [  --with-moc=MOC          use MOC as the absolute path to the QT meta object compiler (moc)])
AC_ARG_WITH(uic,               [  --with-uic=UIC          use UIC as the absolute path to the QT user interface compiler (uic)])
AC_ARG_WITH(opengl-incl,       [  --with-opengl-incl=DIR  OpenGL/Mesa header files are in DIR/GL])
AC_ARG_WITH(opengl-libs,       [  --with-opengl-libs=DIR  OpenGL/Mesa libraries are in DIR/GL])
AC_ARG_WITH(default-path,      [  --with-default-path=DIR use default path DIR/include for contrib header files and DIR/lib for contrib libraries ("../contrib" if unset)])
AC_ARG_WITH(cgal-incl,         [  --with-cgal-incl=DIR    CGAL header files are in DIR])
AC_ARG_WITH(cgal-libs,         [  --with-cgal-libs=DIR    CGAL libraries are in DIR])
AC_ARG_WITH(mesa,              [  --with-mesa             use MESA instead of OpenGL])
AC_ARG_ENABLE(netcdf,          [  --disable-netcdf        disable the use of NetCDF (a library for a data format)])
AC_ARG_WITH(netcdf-incl,       [  --with-netcdf-incl=DIR  specify the location of the NetCDF header file "netcdf.h"])
AC_ARG_WITH(netcdf-libs,       [  --with-netcdf-libs=DIR  specify the location of the NetCDF archive files "libnetcdf.a" and "libnetcdf_c++.a"])
AC_ARG_ENABLE(andims,          [  --disable-andims        disable the use of ANDI/MS (a library for a data format based on NetCDF)])
AC_ARG_WITH(andims-incl,       [  --with-andims-incl=DIR  specify the location of the ANDI/MS header file "ms10.h"])
AC_ARG_WITH(andims-libs,       [  --with-andims-libs=DIR  specify the location of the ANDI/MS archive file "ms10lib.a"])
AC_ARG_WITH(gsl-incl,          [  --with-gsl-incl=DIR     specify the location of the directory containing "gsl/gsl_version.h"])
AC_ARG_WITH(gsl-libs,          [  --with-gsl-libs=DIR     specify the location of the GSL library files "libgsl.so" and "libgslcblas.so"])
AC_ARG_WITH(libsvm-incl,       [  --with-libsvm-incl=DIR  specify the location of the libSVM directory containing "svm.h"])
AC_ARG_WITH(libsvm-libs,       [  --with-libsvm-libs=DIR  specify the location of the libSVM directory containing "svm.o"])
AC_ARG_WITH(xerces-incl,       [  --with-xerces-incl=DIR  specify the location of the XERCES include directory containing "xercesc/util/Base64.hpp"])
AC_ARG_WITH(xerces-libs,       [  --with-xerces-libs=DIR  specify the location of the XERCES lib directory containing "libxerces-c.so"])
AC_ARG_WITH(seqan-incl,        [  --with-seqan-incl=DIR   specify the location of the SEQAN include directory containig the headers])

dnl
dnl find the path to the OpenMS root directory
dnl

if test "x`pwd -W 2>/dev/null`" = "x"; then
  OPENMS_PATH=`cd ..; pwd`
else
  OPENMS_PATH=`cd ..; pwd -W`		# use -W switch for MinGW, as it will give the Windows instead of Posix path
fi


dnl
dnl    interpret switches...
dnl

if test "$prefix" != "NONE" ; then
	PREFIX="${prefix}"
else
	PREFIX="${OPENMS_PATH}"
fi


if test "$with_optimization" = "no" -o "$enable_optimization" = "no"; then
  OPTIMIZATION=false
fi

if test "$enable_debuginfo" = "yes" ; then
  DEBUG_INFO=true
  OPTIMIZATION=false
else
	DEBUG_INFO=false
fi

AC_MSG_CHECKING(STL_debug)
if test "$enable_stl_debug" = "yes" ; then
  GCC_STL_DEBUG=-D_GLIBCXX_DEBUG
	AC_MSG_RESULT([enabled])
else
	AC_MSG_RESULT([disabled])
	GCC_STL_DEBUG=
fi

dnl
dnl   checking for externalmemory-Flag
dnl
AC_MSG_CHECKING(externalmemory)
if test "$enable_externalmemory" = "yes" -o "$with_externalmemory" = "yes"; then
	AC_DEFINE(OPENMS_ENABLE_EXTERNALMEMORY, 1)
	AC_MSG_RESULT(enabled)
else
	AC_MSG_RESULT(disabled (default))
fi

dnl
dnl   checking for openmp-flag
dnl
AC_MSG_CHECKING(openmp)
if test "$enable_openmp" = "yes" -o "$with_openmp" = "yes"; then
	CXXFLAGS="$CXX_FLAGS -fopenmp"
	LIBS="$LIBS -fopenmp"
	AC_MSG_RESULT(enabled)
else
	AC_MSG_RESULT(disabled (default))
fi


if test "$with_pedantic" = "no" -o "$enable_pedantic" = "no"; then
	GCC_PEDANTIC=
else
	GCC_PEDANTIC=-pedantic
fi

if test "$with_default_path" != "" ; then
  DEFAULT_LIB_PATH="${with_default_path}/lib"
  DEFAULT_INCL_PATH="${with_default_path}/include"
else
  DEFAULT_LIB_PATH="${OPENMS_PATH}/contrib/lib"
  DEFAULT_INCL_PATH="${OPENMS_PATH}/contrib/include"
fi

if test "$with_db_plugin" != ""; then
	DB_PLUGIN="${with_db_plugin}"
else
	DB_PLUGIN="QMYSQL"
fi
AC_DEFINE_UNQUOTED(DB_PLUGIN,"${DB_PLUGIN}")

if test "$with_netcdf" = "no" -o "$enable_netcdf" = "no"; then
  NETCDF_SUPPORT=false
else
  NETCDF_SUPPORT=true
  AC_DEFINE(NETCDF_DEF)
fi

if test "$with_netcdf_incl" != "" ; then
  NETCDF_INCPATH="${with_netcdf_incl}"
else
  NETCDF_INCPATH=$DEFAULT_INCL_PATH
fi

if test "$with_netcdf_libs" != "" ; then
  NETCDF_LIBPATH="${with_netcdf_libs}"
else
  NETCDF_LIBPATH=$DEFAULT_LIB_PATH
fi


if test "$with_andims" = "no" -o "$enable_andims" = "no" -o "${NETCDF_SUPPORT}" = false ; then
  ANDIMS_SUPPORT=false
else
  ANDIMS_SUPPORT=true
  AC_DEFINE(ANDIMS_DEF)
fi


if test "$with_andims_incl" != "" ; then
  ANDIMS_INCPATH="${with_andims_incl}"
else
  ANDIMS_INCPATH=$DEFAULT_INCL_PATH
fi

if test "$with_andims_libs" != "" ; then
  ANDIMS_LIBPATH="${with_andims_libs}"
else
  ANDIMS_LIBPATH=$DEFAULT_LIB_PATH
fi

#GSL is always on
GSL_SUPPORT=true
AC_DEFINE(GSL_DEF)

if test "$with_gsl_incl" != "" ; then
  GSL_INCPATH="${with_gsl_incl}"
else
  GSL_INCPATH=$DEFAULT_INCL_PATH
fi

if test "$with_gsl_libs" != "" ; then
  GSL_LIBPATH="${with_gsl_libs}"
else
  GSL_LIBPATH=$DEFAULT_LIB_PATH
fi

#SeqAn is a template library, just needs headers and is always on
AC_DEFINE(SEQAN_DEF)

if test "${with_seqan_incl}" != "" ; then
	SEQAN_INCPATH="${with_seqan_incl}"
else
	SEQAN_INCPATH=$DEFAULT_INCL_PATH
fi

CGAL_SUPPORT=true
AC_DEFINE(CGAL_DEF)

if test "$with_cgal_incl" != "" ; then
  CGAL_INCPATH="${with_cgal_incl}"
else
  CGAL_INCPATH=$DEFAULT_INCL_PATH
fi

if test "$with_cgal_libs" != "" ; then
  CGAL_LIBPATH="${with_cgal_libs}"
else
  CGAL_LIBPATH=$DEFAULT_LIB_PATH
fi

if test "$with_compiler" != "" ; then
	CXX=$with_compiler
fi

if test "$with_cxxflags" != "" ; then
	CXXFLAGS="${CXXFLAGS}`echo $with_cxxflags|tr , \  `"
fi

if test "$with_ldflags" != "" ; then
	LDFLAGS="${LDFLAGS}`echo $with_ldflags|tr ,  \  `"
fi

if test "$with_dynarflags" != "" ; then
	DYNAROPTS="${DYNAROPTS}`echo $with_dynarflags|tr ,  \  `"
fi

if test "$with_arflags" != "" ; then
	AROPTS="${AROPTS}`echo $with_arflags|tr ,  \  `"
fi

#Qt support
QT_INCPATH=""
QT_LIBPATH=""
MOC=moc
UIC=uic
if test "$with_qt" != "" ; then
	QT_INCPATH=$with_qt/include/
	QT_LIBPATH=$with_qt/lib/
	MOC=$with_qt/bin/moc
	UIC=$with_qt/bin/uic
fi
if test "$with_qt_incl" != "" ; then
  QT_INCPATH=$with_qt_incl
fi
if test "$with_qt_libs" != "" ; then
  QT_LIBPATH=$with_qt_libs
fi
if test "$with_moc" != "" ; then
  MOC=$with_moc
fi
if test "$with_uic" != "" ; then
  UIC=$with_uic
fi

#OpenGL/Mesa support
if test "$with_opengl_incl" != "" ; then
  OPENGL_INCPATH=$with_opengl_incl
fi
if test "$with_opengl_libs" != "" ; then
  OPENGL_LIBPATH=$with_opengl_libs
fi
if test "$with_mesa" != "yes" ; then
  GUI_PLATFORM=OpenGL
else
 	GUI_PLATFORM=Mesa
fi

# support for libSVM
LIBSVM_SUPPORT=true
AC_DEFINE(LIBSVM_DEF)
if test "$with_libsvm_incl" != "" ; then
  LIBSVM_INCPATH="${with_libsvm_incl}"
else
  LIBSVM_INCPATH=$DEFAULT_INCL_PATH
fi

if test "$with_libsvm_libs" != "" ; then
  LIBSVM_LIBPATH="${with_libsvm_libs}"
else
  LIBSVM_LIBPATH=$DEFAULT_LIB_PATH
fi

# support for XERCES
XERCES_SUPPORT=true
AC_DEFINE(XERCES_DEF)

if test "$with_xerces_incl" != "" ; then
  XERCES_INCPATH="${with_xerces_incl}"
else
  XERCES_INCPATH=$DEFAULT_INCL_PATH
fi

if test "$with_xerces_libs" != "" ; then
  XERCES_LIBPATH="${with_xerces_libs}"
else
  XERCES_LIBPATH=$DEFAULT_LIB_PATH
fi

dnl
dnl 	global definitions
dnl

dnl
dnl   check for programs used to determine architecture
dnl

AC_PROG_LN_S
CF_MSG_PATH_PROG(UNAME,uname)
CF_MSG_PATH_PROG(CUT,cut)
CF_MSG_PATH_PROG(TR,tr)
CF_MSG_PATH_PROG(AR,ar)
CF_MSG_PATH_PROG(CP,cp,no)
CF_MSG_PATH_PROG(RM,rm,no)
CF_MSG_PATH_PROG(MV,mv,no)
CF_MSG_PATH_PROG(LN,ln,no)
CF_MSG_PATH_PROG(SED,sed,no)
CF_MSG_PATH_PROG(FIND,find,no)
CF_MSG_PATH_PROG(DIFF,diff,no)
CF_MSG_PATH_PROG(TAIL,tail,no)
CF_MSG_PATH_PROG(GREP,grep,no)
CF_MSG_PATH_PROG(FIND,find,no)
CF_MSG_PATH_PROG(EGREP,egrep,no)
CF_MSG_PATH_PROG(CAT,cat,no)
CF_MSG_PATH_PROG(STRIP,strip,no)

AC_SUBST(MV)
AC_SUBST(RM)

dnl
dnl 		Check whether the default shell (/bin/sh) understands "echo -e".
dnl
CF_CHECK_ECHO
CF_CHECK_FIND
AC_SUBST(ECHO)

dnl
dnl			Make sure we have GNU make
dnl
CHECK_GNU_MAKE
if test "${_cv_gnu_make_command+set}" != "set" ; then
	AC_MSG_RESULT(No GNU make found.)
	AC_MSG_RESULT()
	AC_MSG_RESULT(Please install GNU make on your system (www.gnu.org).)
	CF_ERROR(Aborted.)
fi
MAKE=${_cv_gnu_make_command}

AC_SUBST(MAKE)

dnl
dnl   	Determine OS and architecture and all this stuff.
dnl
CF_DETECT_OS


dnl
dnl			Search for the C++ compiler in the (OS-determined) search order.
CF_SEARCH_CXX

dnl
dnl		Determine the object suffix and executable suffix
dnl
SAVE_CXXFLAGS="${CXXFLAGS}"
AC_PROG_CXX
AC_OBJEXT
AC_EXEEXT
CXXFLAGS="${SAVE_CXXFLAGS}"

dnl
dnl  		Try to identify the compiler and determine its version number.
dnl
CXX_IDENTIFIED=false
if test "${CXX_IDENTIFIED}" != true ; then CF_IDENTIFY_INTEL  fi
if test "${CXX_IDENTIFIED}" != true ; then CF_IDENTIFY_GXX  fi
if test "${CXX_IDENTIFIED}" != true ; then CF_IDENTIFY_KAI  fi
if test "${CXX_IDENTIFIED}" != true ; then CF_IDENTIFY_SUN  fi
if test "${CXX_IDENTIFIED}" != true ; then CF_IDENTIFY_SGI  fi
if test "${CXX_IDENTIFIED}" != true ; then CF_IDENTIFY_COMPAQ  fi

dnl
dnl		Default (fall-back) for make depend
dnl
CXX_MAKEDEPEND="makedepend"
MAKEDEP_CXX_OPTS="-f.Dependencies -D__cplusplus --"

dnl	
dnl		Set the compiler-specific options, archiver, linker, etc.
dnl
case "${CXX_NAME}" in
	g++ )				CF_GXX_OPTIONS ;;
	KAI ) 			CF_KAI_OPTIONS ;;
	Intel ) 		CF_INTEL_OPTIONS ;;
	Compaq )		CF_COMPAQ_OPTIONS ;;
	MIPSpro )		CF_MIPSPRO_OPTIONS ;;
	SunCC )			CF_SUNCC_OPTIONS ;;
esac

CF_BUILD_FULL_CXX_NAME


dnl
dnl   checking for OPTIMIZATION-Flag
dnl
CF_CHECK_OPTIMIZATION_FLAG


dnl
dnl   add libs which SeqAn depends on
dnl
CF_SEQAN_DEPLIBS


dnl
dnl		checks for header files
dnl

dnl
dnl		some platforms need to include ieeefp.h
dnl		for the definition of finite
dnl
if test "$OS" = IRIX -o "$OS" = Solaris ; then
	INCLUDE_IEEEFP=true
	OPENMS_INCLUDE_IEEEFP=true
fi

dnl 	check for standard C headers
dnl 	and how to run the preprocessor
AC_HEADER_STDC



dnl
dnl		Check for some C++/STL features
dnl		  -- the numeric_limits class
dnl     -- null template arguments (whether <> is required)
dnl			-- inline functions as template function arguments
dnl 		-- ANSI-compliant iostreams 
dnl			-- sstream vs. strstream
dnl			-- ARM vs. ANSI-style access modifications
dnl
CF_CHECK_NUM_LIMITS
CF_CHECK_TPL_NULL_ARGS
CF_CHECK_INLINE_TPL_ARGS
CF_CHECK_ANSI_IOSTREAM
CF_CHECK_HAS_SSTREAM
CF_CHECK_ARM_ACCESS_MODIFICATION

dnl
dnl		Check for the word lengths of different C/C++ types, define
dnl			the correct replacement types, and determine the endianness of 
dnl			the machine.
dnl
CF_GET_TYPE_SIZES
CF_C_BIGENDIAN

dnl
dnl		Checks for a couple of headers files and whether they are required.
dnl
AC_TYPE_SIZE_T
AC_HEADER_TIME
CF_CHECK_STDINT_H
CF_CHECK_IEEEFP_H
CF_CHECK_VALUES_H
CF_CHECK_LIMITS_H
CF_CHECK_UNISTD_H
CF_CHECK_PROCESS_H
CF_CHECK_TIME_H
CF_CHECK_DIRENT_H
CF_CHECK_DIRECT_H
CF_CHECK_PWD_H
CF_CHECK_IO_H
CF_CHECK_SYS_TIME_H
CF_CHECK_SYS_STAT_H
CF_CHECK_SYS_TIMES_H
CF_CHECK_SYS_TYPES_H
CF_CHECK_SYS_IOCTL_H
CF_CHECK_SYS_SOCKET_H
CF_CHECK_SYS_PARAM_H
CF_CHECK_NETINET_IN_H
CF_CHECK_NETDB_H
CF_CHECK_ARPA_INET_H

CF_CHECK_SYSCONF
CF_CHECK_KILL

dnl
dnl		Check whether vsnprintf is defined and working correctly.
dnl
CF_CHECK_VSNPRINTF


dnl
dnl Check for NetCDF (a library for a data format)
dnl
CF_NETCDF

dnl
dnl Check for ANDI/MS (a library for a data format based on NetCDF)
dnl
CF_ANDIMS

dnl
dnl Check for GSL (a scientific library)
dnl
CF_GSL

dnl
dnl Check for CGAL (a library used for advanced d-dimensional DS & algs)
dnl
CF_CGAL

dnl
dnl Check for SeqAn (a template library for sequence)
dnl
CF_SEQAN

dnl
dnl Check for libSVM
dnl
CF_LIBSVM

dnl
dnl Check for XERCES (a XML parser)
dnl
CF_XERCES

dnl
dnl Check for the visualization/GUI stuff (QT-dependent).
dnl
QT_MIN_VERSION=430
QT_MAX_VERSION=499
QT_RECOMMENDED_VERSION=4.3.3
CF_GUI

dnl
dnl		valgrind support (leak checking for the OpenMS tests)
dnl
CF_VALGRIND

dnl
dnl	Include the default OpenMS include path and uniquify the list.
dnl
OPENMS_INCLUDES="-I${OPENMS_PATH}/include ${OPENMS_INCLUDES}"
CF_UNIQUIFY(${OPENMS_INCLUDES},OPENMS_INCLUDES)

dnl
dnl	Expand BINFMT to include the compiler name
dnl
BINFMT="${BINFMT}-${CXX_NAME}"

LIBS_EXTERNAL=" ${OPENMS_LIBS} ${LIBS} -lm"
OPENMS_LIBS=" -L${PREFIX}/lib/ -lOpenMS ${LIBS_EXTERNAL}"

dnl
dnl   Remove duplicate library entries, but make sure
dnl   that the linker lines don't get messed up for Mac OS
dnl   (-framework XY needs to be kept as it is)
dnl
if test "${OS}" != "Darwin" ; then
  CF_UNIQUIFY(${LIBS_EXTERNAL},LIBS_EXTERNAL)
  CF_UNIQUIFY(${LIBS_EXTERNAL},LIBS_EXTERNAL)
fi



LIBS="${OPENMS_LIBS}"

dnl
dnl   Join all flags
dnl
CXXFLAGS="$CXXFLAGS $DEFINES"

dnl
dnl  We try to remove dependencies from non-OpenMS files
dnl  to speed up compilation.
dnl  This works only for compilers that create one-line dependencies
dnl  (most compilers except for g++). We simply "grep" away all foreign stuff.
dnl  We can only insert ou egrep filter, if the compiler needs
dnl  a redirection of the output (MAKEDEP_CXX_SUFFIX != "")
dnl
dnl  Btw: escaping is FUN! the expression below uses \ to escape the damned shell,
dnl  [] to escape m4, and finally $$ to escape make - the correct expression
dnl  is something like "(${OpenMS_PATH}|^[^/]*$)". Send me an email if got that immediately...
dnl
if test "${HAS_GPLUSPLUS}" = false -a "${MAKEDEP_CXX_SUFFIX}" != "" ; then
	MAKEDEP_CXX_SUFFIX="|egrep \"(${OpenMS_PATH}|^[[^/]]*\$\$)\" ${MAKEDEP_CXX_SUFFIX}"
fi



dnl  
dnl   Neccessary variable substitutions in *.mak
dnl


dnl	PATHS
AC_SUBST(OPENMS_PATH)         dnl  Path to top level directory
AC_SUBST(OPENMS_INCLUDES)     dnl  Includes 


dnl	DEPENDENCIES

MAKEDEP_CXX_OPTS="${MAKEDEP_CXX_OPTS} ${CXXFLAGS}"
if test ${CXX_MAKEDEPEND} = "makedepend" ; then
	MAKEDEP_CXX_OPTS="${MAKEDEP_CXX_OPTS} -- "
fi
	
AC_SUBST(CXX_MAKEDEPEND)        dnl   path to C++ makedepend
AC_SUBST(MAKEDEP_CXX_OPTS)      dnl   options to C++ makedepend
AC_SUBST(MAKEDEP_CXX_SUFFIX)    dnl   optional suffix (for redirection) for C++ makedepend
AC_SUBST(DEPENDENCY_GREPS)      dnl   command line to remove non-OpenMS headers from the dependencies


dnl         COMPILER
AC_SUBST(CXX)               dnl C++ compiler (absolute path)
CC=CXX
AC_SUBST(CC)
AC_SUBST(CXXFLAGS)          dnl C++ compiler flags
AC_SUBST(ADD_CXXFLAGS)      dnl additional C++ compiler falgs (not passed on to makedepend)
AC_SUBST(LIB_CXXFLAGS)			dnl C++ compiler flags used during compilation of the library only
AC_SUBST(NONLIB_CXXFLAGS)		dnl C++ compiler flags used during compilation of applications, tests, examples only
AC_SUBST(HAS_GPLUSPLUS)     dnl true, if C++ compiler is g++
AC_SUBST(STD_CPP_INCLUDES)  dnl standard include paths of C++ compiler
AC_SUBST(MAKEDEP_OPTS)      dnl options for makedep (needed for Linux)
AC_SUBST(OPENMS_TYPENAME)   dnl if the compiler requires the use of typename
                            dnl this variable is set to "typename". "" otherwise.

dnl         LINKER

AC_SUBST(LDFLAGS)           dnl linker flags
AC_SUBST(LIBS)              dnl linker libraries
AC_SUBST(LIBS_EXTERNAL)     dnl linker libraries (external libs only)
AC_SUBST(DEFINES)           dnl defines (part of LDFLAGS)
AC_SUBST(CXXFLAGS_D)        dnl debug flags for compiler 
AC_SUBST(CXXFLAGS_O)        dnl optimization flags for compiler 
AC_SUBST(CPP_MODE_FLAGS)    dnl set to CXXFLAGS_D|CXXFLAGS_O, depending on OpenMS_DEBUG
AC_SUBST(CPP_MODE_FLAGS_NO_OPTIMIZATION)   dnl set to CXXFLAGS_D|"", depending on OpenMS_DEBUG


dnl	ARCHIVER

dnl if nothing has been defined for AR/AROPTS,
dnl use default

if test "${AR}" = "no" ; then
	AR="ar"
fi
if test "${AROPTS}" = "" ; then
	AROPTS="cru"
fi

AC_SUBST(AR)
AC_SUBST(AROPTS)
AC_SUBST(DYNAR)
AC_SUBST(DYNAROPTS)
AC_SUBST(SED_INPLACE_PARAM)


dnl				  MISC

dnl
dnl   ADDITIONAL_COMMENTS is substituted at the end of the
dnl   information text in Makefile and will be printed after
dnl   the library has been built.
dnl   Any hints concerning a special OS/compiler combination 
dnl   might want to print some information here (as is the
dnl   case for the IRIX LD_LBRARY_PATH/LD_LIBRARY64_PATH/LD_LIBRARYN32_PATH)
dnl

if test "${OS}" = MINGW32 ; then
  ADDITIONAL_COMMENTS="PATH"
fi

if test "${ADDITIONAL_COMMENTS}" = "" ; then
	ADDITIONAL_COMMENTS="LD_LIBRARY_PATH"
fi
AC_SUBST(ADDITIONAL_COMMENTS)

dnl
dnl		MINGW32: We want the nasty console window to disappear for TOPPView and INIFileEditor on Windows
dnl

if test "${OS}" = MINGW32 ; then
	LDFLAGS_GUI_MINGW=" -mwindows"
fi
AC_SUBST(LDFLAGS_GUI_MINGW)


dnl    	Define the installation path
dnl     and the binary format string
AC_DEFINE_UNQUOTED(OPENMS_DATA_PATH, "${PREFIX}/share/OpenMS/")
AC_DEFINE_UNQUOTED(OPENMS_BINFMT, "${BINFMT}")

dnl
dnl	Some compiler specific things (standard compliance)
dnl
dnl	Define OPENMS_NULL_TEMPLATE_ARGS
AC_DEFINE_UNQUOTED(OPENMS_NULL_TEMPLATE_ARGS, ${OPENMS_NULL_TEMPLATE_ARGS})

dnl	Define OPENMS_HAS_ANSI_IOSTREAM
if test "$OPENMS_HAS_ANSI_IOSTREAM" = yes ; then
	AC_DEFINE(OPENMS_HAS_ANSI_IOSTREAM)
fi
dnl	Define OPENMS_HAS_SSTREAM
if test "$OPENMS_HAS_SSTREAM" = yes ; then
	AC_DEFINE(OPENMS_HAS_SSTREAM)
fi

dnl	
dnl    Make target "clean" should also remove 
dnl    template repositories. The name of this directory 
dnl    depends on the compiler and is set in TEMPLATE_DIR.
dnl
AC_SUBST(TEMPLATE_DIR)

dnl
dnl    Create files.
dnl

SUBDIRS="CONCEPT SYSTEM MATH DATASTRUCTURES METADATA KERNEL FORMAT FILTERING COMPARISON CHEMISTRY TRANSFORMATIONS ANALYSIS VISUAL APPLICATIONS"
AC_SUBST(SUBDIRS)

BINFMT_PATH=""

AC_PREFIX_DEFAULT(`cd ..; pwd`)
dnl set prefix as AC_PREFIX_DEFAULT gets not expanded in time

if test "x$prefix" = "xNONE" ; then
	if test "x`pwd -W 2>/dev/null`" = "x"; then
	  prefix=`cd ..; pwd`
	else
	  prefix=`cd ..; pwd -W`	# use -W switch for MinGW, as it will give the Windows instead of Posix path
	fi
fi


dnl
dnl Create output files. Move them to temporary location
dnl only, because some of them will require additional
dnl modification (in MULTI_BUILD mode).
dnl
AC_OUTPUT([
	config_defs.mak.tmp:config/config_defs.mak.in 
 	config.mak.tmp:config/config.mak.in 
	common.mak.tmp:config/common.mak.in 
	Makefile.tmp:config/Makefile.in
	Doxyfile.tmp:config/Doxyfile.in
])
CF_MOVE_CONFIG_FILES

dnl
dnl Make sure the dependencies and object lists are (re)built.
dnl
CF_CLEAR_DEP_FILES

dnl
dnl Just give another helpful message and exit.
dnl

AC_MSG_RESULT()
AC_MSG_RESULT()
AC_MSG_RESULT(Your Makefiles have been created. If everything went OK: just type)
AC_MSG_RESULT()
AC_MSG_RESULT(		make lib)
AC_MSG_RESULT()
AC_MSG_RESULT(to build the OpenMS library.)
AC_MSG_RESULT()
