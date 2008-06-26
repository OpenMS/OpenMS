# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------
# $Maintainer: Marc Sturm, Clemens Groepl $
# --------------------------------------------------------------------------

include ../../config_defs.mak

DIRECTORY=TEST/TOPP

# Tests for individual tools
UNCONDITIONAL_TESTS = \
	FeatureLinker_test \
  MapAligner_test \
	ConsensusID_test \
	MapNormalizer_test \
	AdditiveSeries_test \
	IDMerger_test \
	BaselineFilter_test \
	DBImporter_test \
	DBExporter_test \
	Decharger_test \
	DTAExtractor_test \
	FileConverter_test \
	FileFilter_test \
	FileInfo_test \
	FileMerger_test \
	IDFilter_test \
	MascotAdapter_test \
	SequestAdapter_test\
	InspectAdapter_test\
	PepNovoAdapter_test\
	NoiseFilter_test \
	PeakPicker_test \
	Resampler_test \
	SpectraFilter_test \
	InternalCalibration_test\
	TOFCalibration_test\
	FeatureFinder_test\
	FeatureFinderMRM_test\
	TextExporter_test

ifeq ($(LIBSVM_SUPPORT),true)
LIBSVM_TESTS= \
	RTModel_test \
	RTPredict_test
endif

TOPP_TOOL_TESTS = ${UNCONDITIONAL_TESTS} ${LIBSVM_TESTS}

# Tests for processing pipelines using several tools

# Files to remove
CLEANUP = ${RM} -f *.tmp *.tmp_? *.log core.*

# Shortcut to TOPP bin path.
TOPP = ${PREFIX}/bin/


#######################################################################
# DEFAULT TARGET
# Here all TOPP tool tests are executed.
# Afterwards all temporary data files are removed (*.tmp and *.tmp_?)
default: help $(TOPP_TOOL_TESTS)
	@$(ECHO) ""
	@$(ECHO) "-------------------------------------------------------"
	@$(ECHO) ""
	@$(ECHO) "Congratulations!  All TOPP tests executed successfully."
	@$(ECHO) ""
	@ if test -z ${DEBUG}; then $(CLEANUP); fi

clean:
	@$(CLEANUP)
	@${RM} -f *.output //output files are needed for checker so they are not deleted in CLEANUP

#######################################################################
# HELP TARGET
help:
	@$(ECHO) ""
	@$(ECHO) "Main targets are:"
	@$(ECHO) "  - 'make [default]' normal execution"
	@$(ECHO) "  - 'make debug'     debug mode (no deletion of tmp/log files and diff output)"
	@$(ECHO) "  - 'make verbose'   verbose mode (activates tool output to stdout and stderr)"
	@$(ECHO) "  - 'make dv'        verbose and debug mode"
	@$(ECHO) ""
	@$(ECHO) "To execute only single tests with special modes use:"
	@$(ECHO) "  - 'make <target>'           normal execution"
	@$(ECHO) "  - 'make DEBUG=1 <target>'   activates debug mode"
	@$(ECHO) "  - 'make VERBOSE=1 <target>' activates verbose mode"
	@$(ECHO) ""
	@$(ECHO) "-------------------------------------------------------"
	@$(ECHO) ""

#######################################################################
# The FuzzyDiff tool is similar to diff but tolerates slight numeric and whitespace differences.
FUZZYDIFFINI = FuzzyDiff.ini
FUZZYDIFF = $(OPENMS_PATH)/source/APPLICATIONS/UTILS/FuzzyDiff -ini $(FUZZYDIFFINI)

# This target compiles and installs the FuzzyDiff tool if necessary.
make_fuzzydiff:
	@$(ECHO) "The 'FuzzyDiff' tool is similar to 'diff' but tolerates slight numeric and whitespace differences."
	cd $(OPENMS_PATH)/source/APPLICATIONS/UTILS; $(MAKE) FuzzyDiff
	@$(ECHO) ""
	@$(ECHO) "-------------------------------------------------------"
	@$(ECHO) ""

#######################################################################
# DEBUG MODE
# Changes in debug mode are:
# - diff reports the differences (not only whether files are different)
# - temporary files (*.tmp and *.tmp_?) are not deleted after the execution
debug:
	@$(MAKE) DEBUG=1 default

# Diff that allows for slight numeric and whitespace differences
DIFF_NUMERIC = $(FUZZYDIFF) $(1) $(2) $(3)

DIFF_QUIET = $(call DIFF_NUMERIC,-in1 $(PWD)/$1,-in2 $(PWD)/$2,-instance 1)
DIFF_DEBUG = $(call DIFF_NUMERIC,-in1 $(PWD)/$1,-in2 $(PWD)/$2,-instance 2)

# [ WARNING: There must not be tabs at the beginning of the definitions in the next few lines!  CG ]
ifdef DEBUG
DIFF_WITH_OPT = $(DIFF_DEBUG)
else
DIFF_WITH_OPT = $(DIFF_QUIET)
endif

#######################################################################
# VERBOSE MODE
# Changes in verbose mode are:
# - RUN_PROG_OPT will not redirect the console output (stdout and stderr) of the program to a .output file
verbose:
	@$(MAKE) VERBOSE=1 default

# This function runs a TOPP tool with given options.  To invoke the program,
# write '$(call RUN_PROG_OPT,Toolname,options)'.  [ Don't use RUN_PROG_OPT_INTERNAL
# directly, use RUN_PROG_OPT instead.  Depending on the VERBOSE settings,
# RUN_PROG_OPT will append a redirection or not. ]
ifdef DEBUG
RUN_PROG_OPT_INTERNAL = echo "  Executing: $(TOPP)$(1) $(2)"; $(TOPP)$(1) $(2)
else
RUN_PROG_OPT_INTERNAL = $(TOPP)$(1) $(2)
endif

# [ WARNING: There must not be tabs at the beginning of the definitions in the next few lines!  CG ]
ifdef VERBOSE
RUN_PROG_OPT = $(RUN_PROG_OPT_INTERNAL)
else
# extra ( ) are used because the programm call might have its own indirection (e.g. FileInfo)
RUN_PROG_OPT = ( $(RUN_PROG_OPT_INTERNAL) ) > $(1).output 2>&1
endif

#######################################################################
# VERBOSE AND DEBUG MODE
# Changes in verbose mode are:
# - RUN_PROG_OPT will not redirect the console output (stdout and stderr) of the program to a .output file
dv:
	@$(MAKE) VERBOSE=1 DEBUG=1 default

###########################################################################################
#
#   NOW HERE COMES AN...
#
#   I  M  P  O  R  T  A  N  T     N  O  T  E  :
#   
#   1)  READ
#   2)  UNDERSTAND, AND
#   3)  ADOPT
#   
#   THE FOLLOWING USAGE OF THE TEST_FILE_EQUAL MACRO !!!
#
# The TEST_FILE_EQUAL function checks whether two files are equal (up to whitespace).
# Assume the file to be tested normally has the extension .xtn.
# Then, to invoke the test, you would simply write '$(call TEST_FILE_EQUAL,file.xtn)'.
# The name of the other file is 'file.tmp' by default, that is, the basename
# of the second argument with a suffix '.tmp' appended.  Note that .tmp files are
# cleaned up automatically unless DEBUG mode is set.  (How convenient!)
#
# There may be situations where you need to override this default.  Then
# you can use a three-argument version.  '$(call TEST_FILE_EQUAL,file.xtn)' is
# equivalent to '$(call TEST_FILE_EQUAL,file.tmp,file.xtn)'.  Note the order of
# the file arguments, which was chosen to resemble the syntax of the TEST_EQUAL macro
# defined in CONCEPT/ClassTest.h.  (And don't forget to clean up after the test.)
#
# Thank you.
# 
# 
TEST_FILE_EQUAL = \
  if test "$(suffix $2)" == ".tmp"; \
  then \
		$(ECHO) ""; \
    $(ECHO) "WARNING:  Extension .tmp is reserved for generated files.  Maybe the order of arguments in this test is wrong?  Go and fix $@!"; \
    $(ECHO) ""; \
	fi; \
	if $(call DIFF_WITH_OPT,$(if $2,$1,$(basename $1).tmp),$(if $2,$2,$1)); \
	then \
		$(ECHO) "	File $1: OK"; \
	else \
		$(ECHO) "	$@ failed.  File: $1"; \
		exit 1; \
	fi;

###########################################################################################
# Tests for individual TOPP tools.

AdditiveSeries_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "AdditiveSeries:"
	@$(call RUN_PROG_OPT,AdditiveSeries,-ini AdditiveSeries_1_parameters.ini)
	@$(MV) AdditiveSeries_1_gnuplot_tmp.cmd AdditiveSeries_1_gnuplot_cmd.tmp
	@$(MV) AdditiveSeries_1_gnuplot_tmp.dat AdditiveSeries_1_gnuplot_dat.tmp
	@$(MV) AdditiveSeries_1_gnuplot_tmp.err AdditiveSeries_1_gnuplot_err.tmp
	@$(call TEST_FILE_EQUAL,AdditiveSeries_1_output.txt)
	@$(call TEST_FILE_EQUAL,AdditiveSeries_1_gnuplot_cmd.tmp,AdditiveSeries_1_gnuplot.cmd)
	@$(call TEST_FILE_EQUAL,AdditiveSeries_1_gnuplot_dat.tmp,AdditiveSeries_1_gnuplot.dat)
	@$(call TEST_FILE_EQUAL,AdditiveSeries_1_gnuplot_err.tmp,AdditiveSeries_1_gnuplot.err)


IDMerger_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "IDMerger:"
	@$(call RUN_PROG_OPT,IDMerger,-ini IDMerger_1_parameters.ini -out IDMerger_1_output.tmp_)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' IDMerger_1_output.tmp_ > IDMerger_1_output.tmp
	@$(call TEST_FILE_EQUAL,IDMerger_1_output.idXML)


BaselineFilter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "BaselineFilter:"
	@$(call RUN_PROG_OPT,BaselineFilter, -in RawData_input.mzData -out BaselineFilter_output.tmp -resampling 0.2 -struc_elem_length 1.5)
	@$(call TEST_FILE_EQUAL,BaselineFilter_output.mzData)


MapNormalizer_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "MapNormalizer:"
	@$(call RUN_PROG_OPT,MapNormalizer, -in MapNormalizer_input.mzData -out MapNormalizer_output.tmp)
	@$(call TEST_FILE_EQUAL,MapNormalizer_output.mzData)

DBExporter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "DBExporter:"
	@if test -f DBExporter_input.ini; then \
		$(ECHO) "- Test 1: -id"; \
		$(call RUN_PROG_OPT,DBExporter, -out DBExporter_output.tmp -ini DBExporter_input.ini -id 1); \
		($(call TEST_FILE_EQUAL,DBExporter_output.tmp,DBImporter_input.mzData)); \
		$(ECHO) "- Test 2: -query"; \
		$(call RUN_PROG_OPT,DBExporter, -out DBExporter_output2.tmp -ini DBExporter_input.ini -query 'SELECT 1'); \
		($(call TEST_FILE_EQUAL,DBExporter_output2.tmp,DBImporter_input.mzData)); \
	else \
		$(ECHO) "  Cannot run test as there is no INI file 'DBExporter_input.ini'"; \
	fi;

DBImporter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "DBImporter:"
	@if test -f DBImporter_input.ini; then \
		$(ECHO) "- Test 1: -init"; \
		$(call RUN_PROG_OPT,DBImporter, -init -ini DBImporter_input.ini); \
		$(ECHO) "- Test 2: -in"; \
		$(call RUN_PROG_OPT,DBImporter, -in DBImporter_input.mzData -ini DBImporter_input.ini); \
		$(call RUN_PROG_OPT,DBImporter, -in DBImporter_input2.mzData -ini DBImporter_input.ini); \
		$(ECHO) "  The actual test is perfomed by DBExporter!"; \
	else \
		$(ECHO) "  Cannot run test as there is no INI file 'DBImporter_input.ini'"; \
	fi;


Decharger_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "Decharger:"
	@$(call RUN_PROG_OPT,Decharger, -ini Decharger.ini -in Decharger_input.featureXML -out Decharger_output.tmp)
	@$(call TEST_FILE_EQUAL,Decharger_output.featureXML)


DTAExtractor_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "DTAExtractor:"

	@$(ECHO) "- Test 1: rt"
	@$(call RUN_PROG_OPT,DTAExtractor,-in DTAExtractor_1_input.mzData -out DTAExtractor -rt :61)
	@$(MV) DTAExtractor_RT60.dta DTAExtractor.tmp
	@$(call TEST_FILE_EQUAL,DTAExtractor.tmp,DTAExtractor_1_output.dta)

	@$(ECHO) "- Test 2: level"
	@$(call RUN_PROG_OPT,DTAExtractor,-in DTAExtractor_1_input.mzData -out DTAExtractor -level 1)
	@$(MV) DTAExtractor_RT60.dta DTAExtractor.tmp
	@$(call TEST_FILE_EQUAL,DTAExtractor.tmp,DTAExtractor_2_output.dta)

	@$(ECHO) "- Test 3: precursor mz, level"
	@$(call RUN_PROG_OPT,DTAExtractor,-in DTAExtractor_1_input.mzData -out DTAExtractor -level 2 -mz :1000)
	@$(MV) DTAExtractor_RT140_MZ5.dta DTAExtractor.tmp
	@$(call TEST_FILE_EQUAL,DTAExtractor.tmp,DTAExtractor_3_output.dta)

FeatureFinder_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "FeatureFinder:"

	@$(ECHO) "- Test 1: raw data with algorithm 'simple'"
	@$(call RUN_PROG_OPT,FeatureFinder, -ini FeatureFinder_1_parameters.ini -in FeatureFinder_1_input.mzData -out FeatureFinder_1_output.tmp)
	@$(call TEST_FILE_EQUAL,FeatureFinder_1_output.xml)

	@$(ECHO) "- Test 2: raw data with algorithm 'simplest'"
	@$(call RUN_PROG_OPT,FeatureFinder, -ini FeatureFinder_2_parameters.ini -in FeatureFinder_1_input.mzData -out FeatureFinder_2_output.tmp)
	@$(call TEST_FILE_EQUAL,FeatureFinder_2_output.xml)

	@$(ECHO) "- Test 3: raw data with algorithm 'isotope_wavelet'"
	@$(call RUN_PROG_OPT,FeatureFinder, -ini FeatureFinder_3_parameters.ini -in FeatureFinder_1_input.mzData -out FeatureFinder_3_output.tmp)
	@$(call TEST_FILE_EQUAL,FeatureFinder_3_output.xml)

	@$(ECHO) "- Test 4: peak data with algorithm 'picked_peak'"
	@$(call RUN_PROG_OPT,FeatureFinder, -ini FeatureFinder_4_parameters.ini -in FeatureFinder_4_input.mzData -out FeatureFinder_4_output.tmp)
	@$(call TEST_FILE_EQUAL,FeatureFinder_4_output.xml)

FileConverter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "FileConverter:"

	@$(ECHO) "- Test 1: MzData -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_1_input.txt -in_type mzData -out FileConverter.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileConverter.tmp,FileConverter_1_output.dta2d)

	@$(ECHO) "- Test 2: DTA2D -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_2_input.txt -in_type DTA2D -out FileConverter.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileConverter.tmp,FileConverter_2_output.dta2d)

	@$(ECHO) "- Test 3: FeatureFile -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_3_input.txt -in_type featureXML -out FileConverter.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileConverter.tmp,FileConverter_3_output.dta2d)

	@$(ECHO) "- Test 4: MzXML -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_4_input.txt -in_type mzXML -out FileConverter.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileConverter.tmp,FileConverter_4_output.dta2d)

ifeq ($(ANDIMS_SUPPORT),true)
	@$(ECHO) "- Test 5: netCDF -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_5_input.txt -in_type cdf -out FileConverter.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileConverter.tmp,FileConverter_5_output.dta2d)
endif
	@$(ECHO) "- Test 6a: DTA2D -> MzXML -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_6_input.dta2d -out FileConverter.tmp -out_type mzXML)
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter.tmp -in_type mzXML -out FileConverter2.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileConverter2.tmp,FileConverter_6_input.dta2d)

	@$(ECHO) "- Test 6b: DTA2D -> MzData -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_6_input.dta2d -out FileConverter.tmp -out_type mzData)
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter.tmp -in_type mzData -out FileConverter2.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileConverter2.tmp,FileConverter_6_input.dta2d)

	@$(ECHO) "- Test 7: consensusXML -> featureXML"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_7_input.consensusXML -out FileConverter_7_output.tmp -out_type featureXML)
	@$(call TEST_FILE_EQUAL,FileConverter_7_output.tmp,FileConverter_7_output.featureXML)


FileFilter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "FileFilter:"

	@$(ECHO) "- Test 1: RT, m/z test, int (upper bound)"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter_1_out.tmp -rt :30 -mz :1000 -int :20000)
	@$(call TEST_FILE_EQUAL,FileFilter_1_out.tmp,FileFilter_1_output.mzData)

	@$(ECHO) "- Test 2: RT, m/z test, int (lower bound)"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter_2_out.tmp -rt 30: -mz 1000: -int 100:)
	@$(call TEST_FILE_EQUAL,FileFilter_2_out.tmp,FileFilter_2_output.mzData)

	@$(ECHO) "- Test 3: level"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter_3_out.tmp -level 2)
	@$(call TEST_FILE_EQUAL,FileFilter_3_out.tmp,FileFilter_3_output.mzData)

	@$(ECHO) "- Test 4: removing zoom scans"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_4_input.mzData -out FileFilter_4_out.tmp -remove_mode SelectedIonDetection)
	@$(call TEST_FILE_EQUAL,FileFilter_4_out.tmp,FileFilter_4_output.mzData)

	@$(ECHO) "- Test 5: featureXML (upper bound)"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_5_input.featureXML -out FileFilter_5_out.tmp -rt :1000 -mz :480 -int :79000 -charge :3 -q :0.6)
	@$(call TEST_FILE_EQUAL,FileFilter_5_out.tmp,FileFilter_5_out.featureXML)

	@$(ECHO) "- Test 6: featureXML (lower bound)"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_5_input.featureXML -out FileFilter_6_out.tmp -rt 1000: -mz 440: -int 70000: -charge 3: -q 0.51:)
	@$(call TEST_FILE_EQUAL,FileFilter_6_out.tmp,FileFilter_6_out.featureXML)

	@$(ECHO) "- Test 7: Filtering of meta data"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_7_input.mzData -out FileFilter_7_out.tmp -int 7000:)
	@$(call TEST_FILE_EQUAL,FileFilter_7_out.tmp,FileFilter_7_out.mzData)

	@$(ECHO) "- Test 8: consensusXML"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_8_input.consensusXML -out FileFilter_8_out.tmp -rt 600:1400 -mz 700:2300 -int  1100:6000)
	@$(call TEST_FILE_EQUAL,FileFilter_8_out.tmp,FileFilter_8_out.consensusXML)


FileInfo_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "FileInfo:"

	@$(ECHO) "- Test 1: DTA file test"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_1_input.txt -in_type DTA -no_progress > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo.tmp,FileInfo_1_output.txt)

	@$(ECHO) "- Test 2: DTA2D file test"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_2_input.dta2d -no_progress > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo.tmp,FileInfo_2_output.txt)

	@$(ECHO) "- Test 3: feature file test (with meta info and statistics)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_3_input.featureXML -m -s -no_progress > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo.tmp,FileInfo_3_output.txt)

	@$(ECHO) "- Test 4: MzXML file test (with meta info)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_4_input.mzXML -m -no_progress > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo.tmp,FileInfo_4_output.txt)

	@$(ECHO) "- Test 5: MzData file test (with meta info and statistics)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_5_input.txt -in_type mzData -m -s -no_progress > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo.tmp,FileInfo_5_output.txt)

	@$(ECHO) "- Test 6: Meta data arrays test (with details and statistics)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_6_input.mzData -d -s -no_progress > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo.tmp,FileInfo_6_output.txt)

	@$(ECHO) "- Test 7: ConsensusXML file test (with statistics)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_7_input.consensusXML -s -no_progress > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo.tmp,FileInfo_7_output.txt)

ifeq ($(ANDIMS_SUPPORT),true)
	@$(ECHO) "- Test 8: netCDF file test (with meta info)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_8_input.notype -in_type cdf -m -no_progress > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo.tmp,FileInfo_8_output.txt)
endif



FileMerger_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "FileMerger:"

	@$(ECHO) "- Test 1: 3 Files (-rt_file)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_1_ini.txt -out FileMerger.tmp -rt_file)
	@$(call RUN_PROG_OPT,FileConverter,-in FileMerger.tmp -in_type mzData -out FileMerger2.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileMerger2.tmp,FileMerger_1_output.txt)

	@$(ECHO) "- Test 2: 2 Files (-rt_file -in_type)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_2_ini.txt -in_type DTA -out FileMerger.tmp -rt_file)
	@$(call RUN_PROG_OPT,FileConverter,-in FileMerger.tmp -in_type mzData -out FileMerger2.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileMerger2.tmp,FileMerger_2_output.txt)

	@$(ECHO) "- Test 3: 2 Files (-rt_auto)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_3_ini.txt -out FileMerger.tmp -rt_auto)
	@$(call RUN_PROG_OPT,FileConverter,-in FileMerger.tmp -in_type mzData -out FileMerger2.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileMerger2.tmp,FileMerger_3_output.txt)

	@$(ECHO) "- Test 4: 2 Files (several files)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_4_ini.txt -out FileMerger.tmp)
	@$(call RUN_PROG_OPT,FileConverter,-in FileMerger.tmp -in_type mzData -out FileMerger2.tmp -out_type DTA2D)
	@$(call TEST_FILE_EQUAL,FileMerger2.tmp,FileMerger_4_output.txt)

	@$(ECHO) "- Test 5: 3 Files (guess rt from filename, -rt_from_filename, -user_ms_level, -ms_level)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_5_ini.txt -rt_from_filename -user_ms_level -ms_level 2  -out FileMerger5.tmp)
	@$(call TEST_FILE_EQUAL,FileMerger5.tmp,FileMerger_5_output.txt)

IDFilter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "IDFilter:"

	@$(ECHO) "- Test 1: sequence file"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_1_input.idXML -out IDFilter_1_output.tmp_ -sequences_file IDFilter_1_input.fas)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' IDFilter_1_output.tmp_ > IDFilter_1_output.tmp
	@$(call TEST_FILE_EQUAL,IDFilter_1_output.idXML)

	@$(ECHO) "- Test 2: threshold fractions"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_2_input.idXML -out IDFilter_2_output.tmp_ -pep_fraction 1 -prot_fraction 1)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' IDFilter_2_output.tmp_ > IDFilter_2_output.tmp
	@$(call TEST_FILE_EQUAL,IDFilter_2_output.idXML)

	@$(ECHO) "- Test 3: exclusion peptides"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_3_input.idXML -out IDFilter_3_output.tmp_ -exclusion_peptides_file IDFilter_3_2_input.idXML)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' IDFilter_3_output.tmp_ > IDFilter_3_output.tmp
	@$(call TEST_FILE_EQUAL,IDFilter_3_output.idXML)

	@$(ECHO) "- Test 4: retention time"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_4_input.idXML -out IDFilter_4_output.tmp_ -rt_filtering -p_value 0.08)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' IDFilter_4_output.tmp_ > IDFilter_4_output.tmp
	@$(call TEST_FILE_EQUAL,IDFilter_4_output.idXML)

	@$(ECHO) "- Test 5: threshold scores"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_5_input.idXML -out IDFilter_5_output.tmp_ -pep_score 32 -prot_score 25)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' IDFilter_5_output.tmp_ > IDFilter_5_output.tmp
	@$(call TEST_FILE_EQUAL,IDFilter_5_output.idXML)

	@$(ECHO) "- Test 6: best n hits"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_6_input.idXML -out IDFilter_6_output.tmp_ -best_n_peptide_hits 2 -best_n_protein_hits 10)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' IDFilter_6_output.tmp_ > IDFilter_6_output.tmp
	@$(call TEST_FILE_EQUAL,IDFilter_6_output.idXML)

MapAligner_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "MapAligner:"

	@$(ECHO) "- Test 1: pose_clustering feature files"
	@$(call RUN_PROG_OPT,MapAligner,-ini MapAligner_1_parameters.ini)
	@$(call TEST_FILE_EQUAL,MapAligner_1_output1.tmp,MapAligner_1_output1.xml)
	@$(call TEST_FILE_EQUAL,MapAligner_1_output2.tmp,MapAligner_1_output2.xml)
	@$(call TEST_FILE_EQUAL,MapAligner_1_output3.tmp,MapAligner_1_output3.xml)

	@$(ECHO) "- Test 2: pose_clustering peak files"
	@$(call RUN_PROG_OPT,MapAligner,-ini MapAligner_2_parameters.ini)
	@$(call TEST_FILE_EQUAL,MapAligner_2_output1.tmp,MapAligner_2_output1.xml)
	@$(call TEST_FILE_EQUAL,MapAligner_2_output2.tmp,MapAligner_2_output2.xml)
	@$(call TEST_FILE_EQUAL,MapAligner_2_output3.tmp,MapAligner_2_output3.xml)

	@$(ECHO) "- Test 3: spectrum_alignment peak files"
	@$(call RUN_PROG_OPT,MapAligner,-ini MapAligner_3_parameters.ini)
	@$(call TEST_FILE_EQUAL,MapAligner_3_output1.tmp,MapAligner_3_output1.xml)
	@$(call TEST_FILE_EQUAL,MapAligner_3_output2.tmp,MapAligner_3_output2.xml)
	@$(call TEST_FILE_EQUAL,MapAligner_3_output3.tmp,MapAligner_3_output3.xml)

MascotAdapter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "MascotAdapter:"
	@$(call RUN_PROG_OPT,MascotAdapter,-ini MascotAdapter_1_parameters.ini -mascot_in)
	@$(call TEST_FILE_EQUAL,MascotAdapter_1_output.tmp,MascotAdapter_1_output.mascot_in)
	@$(call RUN_PROG_OPT,MascotAdapter,-in MascotAdapter_2_input.mascotXML -out MascotAdapter_2_output.tmp_ -mascot_out)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' MascotAdapter_2_output.tmp_ > MascotAdapter_2_output.tmp
	@$(call TEST_FILE_EQUAL,MascotAdapter_2_output.tmp,MascotAdapter_2_output.idXML)


NoiseFilter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "NoiseFilter:"
	@$(ECHO) "- Test 1: Savitzky-Golay filter"
	@$(call RUN_PROG_OPT,NoiseFilter, -ini NoiseFilter_1_parameters.ini -in RawData_input.mzData -out NoiseFilter_1_output.tmp)
	@$(call TEST_FILE_EQUAL,NoiseFilter_1_output.tmp,NoiseFilter_1_output.mzData)


PeakPicker_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "PeakPicker:"

	@$(ECHO) "-Test1 :PeakPicker without deconvolution"
	@$(call RUN_PROG_OPT,PeakPicker, -ini PeakPicker_parameters.ini -in RawData_input.mzData -out PeakPicker_output.tmp)
	@$(call TEST_FILE_EQUAL,PeakPicker_output.mzData)

	@$(ECHO) ""
	@$(ECHO) "-Test2 :PeakPicker with deconvolution"
	@$(call RUN_PROG_OPT,PeakPicker, -ini PeakPicker_deconv_parameters.ini -in PeakPicker_deconv_input.mzData -out PeakPicker_deconv_output.tmp)
	@$(call TEST_FILE_EQUAL,PeakPicker_deconv_output.mzData)

	@$(ECHO) ""
	@$(ECHO) "-Test3 :PeakPicker with 2D optimization"
	@$(call RUN_PROG_OPT,PeakPicker, -ini PeakPicker_2Dopt_parameters.ini -in RawData_input.mzData -out PeakPicker_2Dopt_output.tmp)
	@$(call TEST_FILE_EQUAL,PeakPicker_2Dopt_output.mzData)


Resampler_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "Resampler:"

	@$(ECHO) "- Test 1: output mzData"
	@$(call RUN_PROG_OPT,Resampler,-in Resampler_1_input.mzData -out Resampler_1_output.tmp -mz 685:697 -cols_mz 20 -rt 1230:1260 -rows_rt 15)
	@$(call TEST_FILE_EQUAL,Resampler_1_output.mzData)

	@$(ECHO) "- Test 2: output png"
	@$(call RUN_PROG_OPT,Resampler,-in Resampler_1_input.mzData -png Resampler_2_output.tmp -cols_mz 20 -rows_rt 15)
	@$(call TEST_FILE_EQUAL,Resampler_2_output.png)


RTModel_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "RTModel:"
	@$(call RUN_PROG_OPT,RTModel,-in RTModel_1_input.idXML -out RTModel_1_output.tmp -ini RTModel_1_parameters.ini)
	@$(call TEST_FILE_EQUAL,RTModel_1_output.model)
	@$(call RUN_PROG_OPT,RTModel,-in_positive RTModel_2_input_positive.idXML -in_negative RTModel_2_input_negative.idXML -out RTModel_2_output.tmp -ini RTModel_2_parameters.ini)
	@$(call TEST_FILE_EQUAL,RTModel_2_output.model)
	@$(call TEST_FILE_EQUAL,RTModel_2_output.tmp_samples,RTModel_2_output.samples)
	@$(call TEST_FILE_EQUAL,RTModel_2_output.tmp_additional_parameters,RTModel_2_output.additional_parameters)
	@$(RM) -rf RTModel_2_output.tmp_additional_parameters RTModel_2_output.tmp_samples

RTPredict_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "RTPredict:"
	@$(call RUN_PROG_OPT,RTPredict,-in RTPredict_1_input.idXML -out RTPredict_1_output.tmp_ -total_gradient_time 3000 -svm_model RTPredict_1_input.model)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' RTPredict_1_output.tmp_ > RTPredict_1_output.tmp
	@$(call TEST_FILE_EQUAL,RTPredict_1_output.idXML)
	@$(call RUN_PROG_OPT,RTPredict,-in RTPredict_2_input.idXML -out_positive RTPredict_2_output_positive.tmp_ -out_negative RTPredict_2_output_negative.tmp_ -svm_model RTPredict_2_input.model)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' RTPredict_2_output_positive.tmp_ > RTPredict_2_output_positive.tmp
	@$(call TEST_FILE_EQUAL,RTPredict_2_output_positive.idXML)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' RTPredict_2_output_negative.tmp_ > RTPredict_2_output_negative.tmp
	@$(call TEST_FILE_EQUAL,RTPredict_2_output_negative.idXML)


InspectAdapter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "InspectAdapter:"

	@$(ECHO) "- Test 1: -list_modifications"
	@$(call RUN_PROG_OPT,InspectAdapter,-ini InspectAdapter_1_parameters.ini -list_modifications > InspectAdapter_1_output.tmp)
	@$(call TEST_FILE_EQUAL,InspectAdapter_1_output.tmp,InspectAdapter_1_output.modification_list.txt)

	@$(ECHO) "- Test 2: -inspect_in with mzXML"
	@$(call RUN_PROG_OPT,InspectAdapter,-ini InspectAdapter_1_parameters.ini -inspect_in -out InspectAdapter_2_output.tmp)
	@$(SED) 's/^spectra,.*$$/spectra,dummy.mzXML/g' InspectAdapter_2_output.tmp > InspectAdapter_2_output.tmp_3
	@$(SED) 's/^db,.*$$/db,dummy.trie/g' InspectAdapter_2_output.tmp_3 > InspectAdapter_2_output.tmp
	@$(call TEST_FILE_EQUAL,InspectAdapter_2_output.tmp,InspectAdapter_2_output.inspect_in)

	@$(ECHO) "- Test 3: -inspect_in with mzData"
	@$(call RUN_PROG_OPT,InspectAdapter,-ini InspectAdapter_1_parameters.ini -inspect_in -in Inspect.mzData -out InspectAdapter_3_output.tmp)
	@$(SED) 's/^spectra,.*$$/spectra,dummy.mzXML/g' InspectAdapter_3_output.tmp > InspectAdapter_3_output.tmp_2
	@$(SED) 's/^db,.*$$/db,dummy.trie/g' InspectAdapter_3_output.tmp_2 > InspectAdapter_3_output.tmp
	@$(call TEST_FILE_EQUAL,InspectAdapter_3_output.tmp,InspectAdapter_2_output.inspect_in)

	@$(ECHO) "- Test 4: -inspect_out"
	@$(call RUN_PROG_OPT,InspectAdapter,-ini InspectAdapter_1_parameters.ini -in InspectAdapter.out -out InspectAdapter_4_output.tmp_ -inspect_out)
	@$(SED) 's/^.*date.*$$/# date removed/g' InspectAdapter_4_output.tmp_ > InspectAdapter_4_output.tmp_2
	@$(SED) 's/^.*xml-stylesheet.*$$//g' InspectAdapter_4_output.tmp_2 > InspectAdapter_4_output.tmp
	@$(call TEST_FILE_EQUAL,InspectAdapter_4_output.tmp,InspectAdapter_4_output.idXML)
	@$(RM) -rf Inspect_FASTAFile_test4.tmp.trie Inspect_FASTAFile_test4.tmp.index InspectAdapter_4_output.tmp_

SequestAdapter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "SequestAdapter:"

	@$(ECHO) "- Test 1: -list_modifications"
	@$(call RUN_PROG_OPT,SequestAdapter,-ini SequestAdapter_1_parameters.ini -list_modifications > SequestAdapter_1_output.tmp)
	@$(call TEST_FILE_EQUAL,SequestAdapter_1_output.tmp,SequestAdapter_1_output.modification_list.txt)

	@$(RM) -rf ./tmp/Sequest*.dta_*  ./tmp/sequest_run.bat ./tmp/temp.sequest.in
	@$(ECHO) "- Test 2: -sequest_in with mzXML"
	@$(call RUN_PROG_OPT,SequestAdapter,-ini SequestAdapter_1_parameters.ini -in Sequest.mzXML -out SequestAdapter_2_output.tmp -sequest_in)
	@$(call TEST_FILE_EQUAL,SequestAdapter_2_output.tmp,SequestAdapter_2_output.sequest_in)
	@$(RM) -rf ./tmp/Sequest*.dta_*

	@$(RM) -rf ./tmp/Sequest*.dta_*  ./tmp/sequest_run.bat ./tmp/temp.sequest.in
	@$(ECHO) "- Test 3: -sequest_in with mzData"
	@$(call RUN_PROG_OPT,SequestAdapter,-ini SequestAdapter_1_parameters.ini -in Sequest.mzData -out SequestAdapter_3_output.tmp -sequest_in -in Sequest.mzData)
	@$(call TEST_FILE_EQUAL,SequestAdapter_3_output.tmp,SequestAdapter_2_output.sequest_in)
	@$(RM) -rf ./tmp/Sequest*.dta_*

	@$(ECHO) "- Test 4: -sequest_out"
	@$(call RUN_PROG_OPT,SequestAdapter,-ini SequestAdapter_2_parameters.ini -in tmp -out SequestAdapter_4_output.tmp_ -sequest_out)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' SequestAdapter_4_output.tmp_ > SequestAdapter_4_output.tmp
	@$(call TEST_FILE_EQUAL,SequestAdapter_4_output.tmp,SequestAdapter_4_output.idXML)
	@$(RM) -rf ./tmp/sequest_run.bat ./tmp/temp.sequest.in

PepNovoAdapter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "PepNovoAdapter:"

	@$(ECHO) "- Test 1: -list_modifications"
	@$(call RUN_PROG_OPT,PepNovoAdapter,-ini PepNovoAdapter_1_parameters.ini -list_modifications > PepNovoAdapter_1_output.tmp)
	@$(call TEST_FILE_EQUAL,PepNovoAdapter_1_output.tmp,PepNovoAdapter_1_output.modification_list.txt)

	@$(ECHO) "- Test 2: -list_models"
	@$(call RUN_PROG_OPT,PepNovoAdapter,-ini PepNovoAdapter_1_parameters.ini -list_models > PepNovoAdapter_2_output.tmp)
	@$(call TEST_FILE_EQUAL,PepNovoAdapter_2_output.tmp,PepNovoAdapter_2_output.model_list.txt)

	@$(ECHO) "- Test 3: -pepnovo_in with mzXML"
	@$(call RUN_PROG_OPT,PepNovoAdapter,-ini PepNovoAdapter_1_parameters.ini -in PepNovo.mzXML -pepnovo_in -out PepNovoAdapter_3_output.tmp)
	@$(call TEST_FILE_EQUAL,./tmp/PepNovo_PTMs_.txt,./tmp/PepNovo_PTMs.txt)
	@$(RM) -rf ./tmp/PepNovo*.dta

	@$(ECHO) "- Test 4: -pepnovo_in with mzData"
	@$(call RUN_PROG_OPT,PepNovoAdapter,-ini PepNovoAdapter_1_parameters.ini -in PepNovo.mzData -pepnovo_in -out PepNovoAdapter_4_output.tmp)
	@$(call TEST_FILE_EQUAL,./tmp/PepNovo_PTMs_.txt,./tmp/PepNovo_PTMs.txt)
	@$(RM) -rf ./tmp/PepNovo*.dta ./tmp/PepNovo_PTMs.txt ./tmp/dta_list.txt

	@$(ECHO) "- Test 5: -pepnovo_out"
	@$(call RUN_PROG_OPT,PepNovoAdapter,-ini PepNovoAdapter_5_parameters.ini -in PepNovoAdapter_5_output.pepnovo_out -out PepNovoAdapter_5_output.tmp_ -pepnovo_out)
	@$(SED) 's/^.*xml-stylesheet.*$$//g' PepNovoAdapter_5_output.tmp_ > PepNovoAdapter_5_output.tmp_2
	@$(SED) 's/^.*date_group_1.*$$/# date group removed/g' PepNovoAdapter_5_output.tmp_2 > PepNovoAdapter_5_output.tmp
	@$(call TEST_FILE_EQUAL,PepNovoAdapter_5_output.tmp,PepNovoAdapter_5_output.idXML)

SpectraFilter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "SpectraFilter:"

	@$(ECHO) "- Test 1: -filter SqrtMower"
	@$(call RUN_PROG_OPT,SpectraFilter,-in SpectraFilter_1_input.mzData -out SpectraFilter_1_output.mzData.tmp -type SqrtMower)
	@$(call TEST_FILE_EQUAL,SpectraFilter_1_output.mzData.tmp,SpectraFilter_1_output.mzData)

	@$(ECHO) "- Test 2: -filter WindowMower"
	@$(call RUN_PROG_OPT,SpectraFilter,-in SpectraFilter_2_input.mzData -out SpectraFilter_2_output.mzData.tmp -type WindowMower)
	@$(call TEST_FILE_EQUAL,SpectraFilter_2_output.mzData.tmp,SpectraFilter_2_output.mzData)

	@$(ECHO) "- Test 3: -filter WindowMower -ini SpectraFilter_3_parameters.ini"
	@$(call RUN_PROG_OPT,SpectraFilter,-in SpectraFilter_3_input.mzData -out SpectraFilter_3_output.mzData.tmp -ini  SpectraFilter_3_parameters.ini -type WindowMower)
	@$(call TEST_FILE_EQUAL,SpectraFilter_3_output.mzData.tmp,SpectraFilter_3_output.mzData)

ConsensusID_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "ConsensusID:"

	@$(ECHO) "- Test 1: Map consensus"
	@$(call RUN_PROG_OPT,ConsensusID,-ini ConsensusID_1_parameters.ini -out ConsensusID_1_output.tmp_)
	@$(SED) 's/^.*IdentificationRun date.*$$/#Run removed because of date/g' ConsensusID_1_output.tmp_ > ConsensusID_1_output.tmp_2
	@$(SED) 's/^.*xml-stylesheet.*$$//g' ConsensusID_1_output.tmp_2 > ConsensusID_1_output.tmp
	@$(call TEST_FILE_EQUAL,ConsensusID_1_output.tmp,ConsensusID_1_output.idXML)

	@$(ECHO) "- Test 2: Feature consensus"
	@$(call RUN_PROG_OPT,ConsensusID,-ini ConsensusID_2_parameters.ini -out ConsensusID_2_output.tmp_ -features ConsensusID_2_features.featureXML)
	@$(SED) 's/^.*IdentificationRun date.*$$/#Run removed because of date/g' ConsensusID_2_output.tmp_ > ConsensusID_2_output.tmp_2
	@$(SED) 's/^.*xml-stylesheet.*$$//g' ConsensusID_2_output.tmp_2 > ConsensusID_2_output.tmp
	@$(call TEST_FILE_EQUAL,ConsensusID_2_output.tmp,ConsensusID_2_output.idXML)

InternalCalibration_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "InternalCalibration:"
	@$(ECHO) "- Test1: calibration of raw data"
	@$(call RUN_PROG_OPT,InternalCalibration,-in InternalCalibration_1_input.mzData -out InternalCalibration_1_output.tmp -ref_masses InternalCalibration_ref_masses -ini InternalCalibration_parameters.ini)
	@$(call TEST_FILE_EQUAL,InternalCalibration_1_output.mzData)
	@$(ECHO) "- Test2: calibration of peak data"
	@$(call RUN_PROG_OPT,InternalCalibration,-in InternalCalibration_2_input.mzData -out InternalCalibration_2_output.tmp -ref_masses InternalCalibration_ref_masses -peak_data)
	@$(call TEST_FILE_EQUAL,InternalCalibration_2_output.mzData)

TOFCalibration_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "TOFCalibration_test:"
	@$(ECHO) "- Test1: calibration of raw data"
	@$(call RUN_PROG_OPT,TOFCalibration,-in TOFCalibration_1_input.mzData -out TOFCalibration_1_output.tmp -ref_masses TOFCalibration_ref_masses -ini TOFCalibration_parameters.ini -tof_const TOFCalibration_const -ext_calibrants TOFCalibration_1_calibrants.mzData)
	@$(call TEST_FILE_EQUAL,TOFCalibration_1_output.mzData)
	@$(ECHO) "- Test2: calibration of peak data"
	@$(call RUN_PROG_OPT,TOFCalibration,-in TOFCalibration_2_input.mzData -out TOFCalibration_2_output.tmp -ref_masses TOFCalibration_ref_masses  -ini TOFCalibration_parameters.ini -tof_const TOFCalibration_const -ext_calibrants TOFCalibration_2_calibrants.mzData -peak_data)
	@$(call TEST_FILE_EQUAL,TOFCalibration_2_output.mzData)

TextExporter_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "TextExporter:"

	@$(ECHO) "- Test 1: FeatureXML test"
	@$(call RUN_PROG_OPT,TextExporter,-in TextExporter_1_input.featureXML -no_progress -out TextExporter.tmp)
	@$(call TEST_FILE_EQUAL,TextExporter.tmp,TextExporter_1_output.txt)

	@$(ECHO) "- Test 2: ConsensusXML test"
	@$(call RUN_PROG_OPT,TextExporter,-in TextExporter_2_input.consensusXML -no_progress -out TextExporter.tmp)
	@$(call TEST_FILE_EQUAL,TextExporter.tmp,TextExporter_2_output.txt)

	@$(ECHO) "- Test 3: IdXML test"
	@$(call RUN_PROG_OPT,TextExporter,-in TextExporter_3_input.idXML -no_progress -out TextExporter3.tmp)
	@$(call TEST_FILE_EQUAL,TextExporter3.tmp,TextExporter_3_output.txt)
	@$(call RUN_PROG_OPT,TextExporter,-in TextExporter_3_input.idXML -no_progress -out TextExporter3.tmp -proteins_only)
	@$(call TEST_FILE_EQUAL,TextExporter3.tmp,TextExporter_3_output_proteins.txt)
	@$(call RUN_PROG_OPT,TextExporter,-in TextExporter_3_input.idXML -no_progress -out TextExporter3.tmp -peptides_only)
	@$(call TEST_FILE_EQUAL,TextExporter3.tmp,TextExporter_3_output_peptides.txt)


FeatureFinderMRM_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "FeatureFinderMRM:"

	@$(ECHO) "- Test 1: MRM quantitation without debug output"
	@$(call RUN_PROG_OPT,FeatureFinderMRM,-in FeatureFinderMRM_1_input.mzData -ini FeatureFinderMRM_1.ini -out FeatureFinderMRM.tmp)
	@$(call TEST_FILE_EQUAL,FeatureFinderMRM.tmp,FeatureFinderMRM_1_output.featureXML)

	@$(ECHO) "- Test 2: MRM quantitation with debug output"
	@$(call RUN_PROG_OPT,FeatureFinderMRM,-in FeatureFinderMRM_1_input.mzData -ini FeatureFinderMRM_1.ini -d -out FeatureFinderMRM.tmp)
	@$(call TEST_FILE_EQUAL,FeatureFinderMRM.tmp,FeatureFinderMRM_1_output.featureXML)
	@$(call TEST_FILE_EQUAL,precursor_461.743,FeatureFinderMRM_2_output.txt)
	@$(RM) -rf precursor_461.743

FeatureLinker_test: make_fuzzydiff
	@$(ECHO) ""
	@$(ECHO) "FeatureLinker:"

	@$(ECHO) "- Test 1: labeled"
	@$(call RUN_PROG_OPT,FeatureLinker,-ini FeatureLinker_1_parameters.ini -in FeatureLinker_1_input.featureXML -out FeatureLinker_1_output.tmp)
	@$(call TEST_FILE_EQUAL,FeatureLinker_1_output.tmp,FeatureLinker_1_output.consensusXML)

	@$(ECHO) "- Test 2: unlabeled (big dataset)"
	@$(call RUN_PROG_OPT,FeatureLinker,-ini FeatureLinker_2_parameters.ini)
	@$(call TEST_FILE_EQUAL,FeatureLinker_2_output.tmp,FeatureLinker_2_output.consensusXML)

	@$(ECHO) "- Test 3: unlabeled (small dataset)"
	@$(call RUN_PROG_OPT,FeatureLinker,-ini FeatureLinker_3_parameters.ini)
	@$(call TEST_FILE_EQUAL,FeatureLinker_3_output.tmp,FeatureLinker_3_output.consensusXML)

	@$(ECHO) "- Test 4: unlabeled (small dataset)"
	@$(call RUN_PROG_OPT,FeatureLinker,-ini FeatureLinker_4_parameters.ini)
	@$(call TEST_FILE_EQUAL,FeatureLinker_4_output.tmp,FeatureLinker_4_output.consensusXML)




