# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

# Tests for individual tools.
# (in alphabetical order)
TOPP_TOOL_TESTS = \
	AdditiveSeries_test \
	AnalysisXMLMerger_test \
	BaselineFilter_test \
	DBExporter_test \
	DBImporter_test \
	DTAExtractor_test \
	Dewarper_test \
	FeatureFinder_test \
	FeaturePairSplitter_test \
	FileConverter_test \
	FileFilter_test \
	FileInfo_test \
	FileMerger_test \
	IDFilter_test \
	InspectAdapter_test \
	LabeledMatcher_test \
	MapAlignment_test \
	MapMatcher_test \
	MapStatistics_test \
	MascotAdapter_test \
	NoiseFilter_test \
	PeakPicker_test \
	RTModel_test \
	RTPredict_test \
	SequestAdapter_test \
	SpectraFilter_test \
	TOPPView_test \
	UnlabeledMatcher_test \


# Tests for processing pipelines using several tools.
# (in alphabetical order)
TOPP_PIPELINE_TESTS = \
	MyoQuant_test \

#files to remove
CLEANUP = ${RM} -f *.tmp  *.rounded_tmp *.log

# Shortcut to TOPP bin path.
TOPP = ${OPENMS_PATH}/bin/

#######################################################################
# DEFAULT TARGET
# Here all TOPP tool tests are executed.
# Afterwards all temporary data files are removed (*.tmp)
default:  help RoundCat $(TOPP_TOOL_TESTS) $(TOPP_PIPELINE_TESTS)
	@$(ECHO) ""
	@$(ECHO) "-------------------------------------------------------"
	@$(ECHO) ""
	@$(ECHO) "Congratulations!  All TOPP tests executed successfully."
	@$(ECHO) ""
	@ if test -z ${DEBUG}; then $(CLEANUP); fi

clean:
	$(CLEANUP)

#######################################################################
# HELP TARGET
help:
	@$(ECHO) ""
	@$(ECHO) "Main targets are:"
	@$(ECHO) "  - 'make [default]' normal execution"
	@$(ECHO) "  - 'make debug'     debug mode (no deletion of tmp/log files and diff output)"
	@$(ECHO) "  - 'make verbose'   verbose mode (activates tool output to stdout and stderr)"
	@$(ECHO) ""
	@$(ECHO) "To execute only single tests with special modes use:"
	@$(ECHO) "  - 'make DEBUG=1 <target>' activates debug mode"
	@$(ECHO) "  - 'make VERBOSE=1 <target>' activates verbose mode"
	@$(ECHO) ""
	@$(ECHO) "-------------------------------------------------------"
	@$(ECHO) ""

#######################################################################
# DEBUG MODE
# Changes in debug mode are:
# - diff reports the differences (not only whether files are different)
# - temporary files (*.tmp) are not deleted after the execution
debug: DEBUG=1 default

# TODO: read from environment _and_ provide a default ????
PRECISION = 3

ROUNDCAT = RoundCat

# Diff that rounds numbers before comparing the files.
DIFF_ROUNDCAT = \
	$(ROUNDCAT) $(1) $(1).rounded_tmp $(PRECISION) && \
	$(ROUNDCAT) $(2) $(2).rounded_tmp $(PRECISION) && \
	$(DIFF) $(1).rounded_tmp $(2).rounded_tmp $(3)

# Diff that ignores whitespaces and empty lines
DIFF_DEBUG = $(call DIFF_ROUNDCAT,$1,$2,-bB --initial-tab -u)
DIFF_QUIET = $(call DIFF_ROUNDCAT,$1,$2,-bB -t --brief)

# [ WARNING: There must not be tabs at the beginning of the definitions in the next few lines!  CG ]
ifdef DEBUG
DIFF_WITH_OPT = $(DIFF_DEBUG)
else
DIFF_WITH_OPT = $(DIFF_QUIET)
endif

#######################################################################
# VERBOSE MODE
# Changes in verbose mode are:
# - RUN_PROG_OPT will not redirect the console output (stdout and stderr) of the program to /dev/null
verbose: VERBOSE=1 default

# This function runs a TOPP tool with given options.  To invoke the program,
# write '$(call RUN_PROG_OPT,Toolname,options)'.  [ Don't use RUN_PROG_OPT_INTERNAL
# directly, use RUN_PROG_OPT instead.  Depending on the VERBOSE settings,
# RUN_PROG_OPT will append a redirection or not. ]
RUN_PROG_OPT_INTERNAL = $(TOPP)$(1) $(2)

# [ WARNING: There must not be tabs at the beginning of the definitions in the next few lines!  CG ]
ifdef VERBOSE
RUN_PROG_OPT = $(RUN_PROG_OPT_INTERNAL)
else
# extra ( ) are used because the programm call might have its own indirection (e.g. FileInfo)
RUN_PROG_OPT = ( $(RUN_PROG_OPT_INTERNAL) ) > /dev/null 2>&1
endif

###########################################################################################
# This function checks whether two files are equal (up to whitespace).  To
# invoke the test, write '$(call TEST_FILE_EQUAL,file.xtn)'.  The name of the
# other file is 'file.tmp' by default, i.e. the basename of the second
# argument with a suffix '.tmp' appended.  To override the default, you can
# specify a third argument.  Thus '$(call TEST_FILE_EQUAL,file.xtn)' is
# equivalent to '$(call TEST_FILE_EQUAL,file.xtn,file.tmp)'.
TEST_FILE_EQUAL =  \
	if $(call DIFF_WITH_OPT,$(if $2,$2,$(basename $1).tmp),$1); \
	then \
		$(ECHO) "	File $1: OK"; \
	else \
		$(ECHO) "	$@ failed.  File: $1"; \
		exit 1; \
	fi;

# A simple minded cat-like application that reads the input file and writes it
# to the output file, but with all numbers rounded to precision digits.
$(ROUNDCAT): $(ROUNDCAT).C

###########################################################################################
# Tests for individual TOPP tools.

AdditiveSeries_test:
	@$(ECHO) ""
	@$(ECHO) "AdditiveSeries:"
	@$(call RUN_PROG_OPT,AdditiveSeries,-ini AdditiveSeries_1_parameters.ini)
	@$(MV) AdditiveSeries_1_gnuplot_tmp.cmd AdditiveSeries_1_gnuplot_cmd.tmp
	@$(MV) AdditiveSeries_1_gnuplot_tmp.dat AdditiveSeries_1_gnuplot_dat.tmp
	@$(MV) AdditiveSeries_1_gnuplot_tmp.err AdditiveSeries_1_gnuplot_err.tmp
	@$(call TEST_FILE_EQUAL,AdditiveSeries_1_output.txt)
	@$(call TEST_FILE_EQUAL,AdditiveSeries_1_gnuplot.cmd,AdditiveSeries_1_gnuplot_cmd.tmp)
	@$(call TEST_FILE_EQUAL,AdditiveSeries_1_gnuplot.dat,AdditiveSeries_1_gnuplot_dat.tmp)
	@$(call TEST_FILE_EQUAL,AdditiveSeries_1_gnuplot.err,AdditiveSeries_1_gnuplot_err.tmp)


AnalysisXMLMerger_test:
	@$(ECHO) ""
	@$(ECHO) "AnalysisXMLMerger:"
	@$(ECHO) "	... todo ..."


BaselineFilter_test:
	@$(ECHO) ""
	@$(ECHO) "BaselineFilter:"
	@$(ECHO) "	... todo ..."


DBExporter_test:
	@$(ECHO) ""
	@$(ECHO) "DBExporter:"
	@$(ECHO) "	... todo ..."


DBImporter_test:
	@$(ECHO) ""
	@$(ECHO) "DBImporter:"
	@$(ECHO) "	... todo ..."


DTAExtractor_test:
	@$(ECHO) ""
	@$(ECHO) "DTAExtractor:"

	@$(ECHO) "- Test 1: rt"
	@$(call RUN_PROG_OPT,DTAExtractor,-in DTAExtractor_1_input.mzData -out DTAExtractor -rt :61)
	@$(MV) DTAExtractor_RT60.dta DTAExtractor.tmp
	@$(call TEST_FILE_EQUAL,DTAExtractor_1_output.dta,DTAExtractor.tmp)

	@$(ECHO) "- Test 2: level"
	@$(call RUN_PROG_OPT,DTAExtractor,-in DTAExtractor_1_input.mzData -out DTAExtractor -level 1)
	@$(MV) DTAExtractor_RT60.dta DTAExtractor.tmp
	@$(call TEST_FILE_EQUAL,DTAExtractor_2_output.dta,DTAExtractor.tmp)

	@$(ECHO) "- Test 3: precursor mz, level"
	@$(call RUN_PROG_OPT,DTAExtractor,-in DTAExtractor_1_input.mzData -out DTAExtractor -level 2 -mz :1000)
	@$(MV) DTAExtractor_RT140_MZ5.dta DTAExtractor.tmp
	@$(call TEST_FILE_EQUAL,DTAExtractor_3_output.dta,DTAExtractor.tmp)


Dewarper_test:
	@$(ECHO) ""
	@$(ECHO) "Dewarper:"
	@$(ECHO) "	... todo ..."


FeatureFinder_test:
	@$(ECHO) ""
	@$(ECHO) "FeatureFinder:"

	@$(ECHO) "- Test 1:"
	@$(call RUN_PROG_OPT,FeatureFinder,\
	-ini FeatureFinder_1_parameters.ini \
	-in FeatureFinder_1_input.mzData \
	-out FeatureFinder_1_output.tmp)
	@$(call TEST_FILE_EQUAL,FeatureFinder_1_output.xml)


FeaturePairSplitter_test:
	@$(ECHO) ""
	@$(ECHO) "FeaturePairSplitter:"
	@$(call RUN_PROG_OPT,FeaturePairSplitter,\
	-in FeaturePairSplitter_1_pairs.xml \
	-out1 FeaturePairSplitter_1_out1.tmp \
	-out2 FeaturePairSplitter_1_out2.tmp \
	-qual FeaturePairSplitter_1_qual.tmp \
	-dump FeaturePairSplitter_1_dump.tmp)
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_out1.xml)
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_out2.xml)
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_qual.wsv)
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_dump.wsv)

FileConverter_test:
	@$(ECHO) ""
	@$(ECHO) "FileConverter:"
	@$(ECHO) "	... todo ..."


FileFilter_test:
	@$(ECHO) ""
	@$(ECHO) "FileFilter:"

	@$(ECHO) "- Test 1: RT, m/z test, int (upper bound)"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter.tmp -rt :30 -mz :1000 -int :20000)
	@$(call TEST_FILE_EQUAL,FileFilter_1_output.mzData,FileFilter.tmp)

	@$(ECHO) "- Test 2: RT, m/z test, int (lower bound)"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter.tmp -rt 30: -mz 1000: -int 100:)
	@$(call TEST_FILE_EQUAL,FileFilter_2_output.mzData,FileFilter.tmp)

	@$(ECHO) "- Test 3: level"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter.tmp -level 2)
	@$(call TEST_FILE_EQUAL,FileFilter_3_output.mzData,FileFilter.tmp)

	@$(ECHO) "- Test 4: remove_zoom"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_4_input.mzData -out FileFilter.tmp -remove_zoom)
	@$(call TEST_FILE_EQUAL,FileFilter_4_output.mzData,FileFilter.tmp)


FileInfo_test:
	@$(ECHO) ""
	@$(ECHO) "FileInfo:"

	@$(ECHO) "- Test 1: DTA file test"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_1_input.txt -in_type DTa > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo_1_output.txt,FileInfo.tmp)

	@$(ECHO) "- Test 2: DTA2D file test"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_2_input.dta2d > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo_2_output.txt,FileInfo.tmp)

	@$(ECHO) "- Test 3: feature file test (with meta info)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_3_input.feat -m > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo_3_output.txt,FileInfo.tmp)

	@$(ECHO) "- Test 4: MzXML file test (with meta info)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_4_input.mzXML -m > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo_4_output.txt,FileInfo.tmp)

	@$(ECHO) "- Test 5: MzData file test (with meta info)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_5_input.txt -in_type MzDaTa -m > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo_5_output.txt,FileInfo.tmp)

FileMerger_test:
	@$(ECHO) ""
	@$(ECHO) "FileMerger:"
	@$(ECHO) "	... todo ..."


IDFilter_test:
	@$(ECHO) ""
	@$(ECHO) "IDFilter:"
	@$(ECHO) "	... todo ..."


InspectAdapter_test:
	@$(ECHO) ""
	@$(ECHO) "InspectAdapter:"
	@$(ECHO) "	... todo ..."


LabeledMatcher_test:
	@$(ECHO) ""
	@$(ECHO) "LabeledMatcher:"
	@$(ECHO) "	... todo ..."


MapAlignment_test:
	@$(ECHO) ""
	@$(ECHO) "MapAlignment:"
	@$(ECHO) "	... todo ..."


MapMatcher_test:
	@$(ECHO) ""
	@$(ECHO) "MapMatcher:"
	@$(ECHO) "	... todo ..."


MapStatistics_test:
	@$(ECHO) ""
	@$(ECHO) "MapStatistics:"
	@$(ECHO) "	... todo ..."


MascotAdapter_test:
	@$(ECHO) ""
	@$(ECHO) "MascotAdapter:"
	@$(ECHO) "	... todo ..."


NoiseFilter_test:
	@$(ECHO) ""
	@$(ECHO) "NoiseFilter:"
	@$(ECHO) "	... todo ..."


PeakPicker_test:
	@$(ECHO) ""
	@$(ECHO) "PeakPicker:"
	@$(ECHO) "	... todo ..."


RTModel_test:
	@$(ECHO) ""
	@$(ECHO) "RTModel:"
	@$(ECHO) "	... todo ..."


RTPredict_test:
	@$(ECHO) ""
	@$(ECHO) "RTPredict:"
	@$(ECHO) "	... todo ..."


SequestAdapter_test:
	@$(ECHO) ""
	@$(ECHO) "SequestAdapter:"
	@$(ECHO) "	... todo ..."


SpectraFilter_test:
	@$(ECHO) ""
	@$(ECHO) "SpectraFilter:"
	@$(ECHO) "	... todo ..."


TOPPView_test:
	@$(ECHO) ""
	@$(ECHO) "TOPPView:"
	@$(ECHO) "	... todo ..."


UnlabeledMatcher_test:
	@$(ECHO) ""
	@$(ECHO) "UnlabeledMatcher:"
	@$(call RUN_PROG_OPT,UnlabeledMatcher,-ini UnlabeledMatcher_1_parameters.ini)
	@$(call TEST_FILE_EQUAL,UnlabeledMatcher_1_pairs.xml)
	@$(call TEST_FILE_EQUAL,UnlabeledMatcher_1_grid.xml)


###########################################################################################
## Tests for processing pipelines.

PRECISION=2

MyoQuant_test:
	@$(ECHO) ""
	@$(ECHO) "MyoQuant_test:"

	@$(ECHO) "- Test 1: Additive Series using intensity ratio of tryptic peptides T11hu / T10ho, no dewarping"
	@ for i in 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37; do \
		$(ECHO) "	Finding features in LC-MS map $$i"; \
		$(call RUN_PROG_OPT,FeatureFinder,-ini MyoQuant_1_parameters.ini -in MyoQuant_1_input$${i}.mzData -out MyoQuant_1_feat$${i}_xml.tmp); \
	done;
	@$(ECHO) "	AdditiveSeries, linear regression";
	@$(call RUN_PROG_OPT,AdditiveSeries,-ini MyoQuant_1_parameters.ini)
	@ $(call TEST_FILE_EQUAL,MyoQuant_1_output.xml)

	@$(ECHO) ""
	@$(ECHO) "- Test 2: Additive Series using intensity ratio of tryptic peptides T11hu / T10ho, with dewarping"
	@$(ECHO) "	... todo ..."


