# -*- Mode: makefile; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
UNCONDITIONAL_TESTS = \
	AdditiveSeries_test \
	AnalysisXMLMerger_test \
	BaselineFilter_test \
	DBImporter_test \
	DBExporter_test \
	DTAExtractor_test \
	MapDewarper_test \
	FeaturePairSplitter_test \
	FileConverter_test \
	FileFilter_test \
	FileInfo_test \
	FileMerger_test \
	IDFilter_test \
	InspectAdapter_test \
	LabeledMatcher_test \
	MapMatcher_test \
	MapStatistics_test \
	MascotAdapter_test \
	NoiseFilter_test \
	SequestAdapter_test\
	SpectraFilter_test \
	FeatureFinder_test \
	PeakPicker_test

#	Resampler_test \

ifeq ($(CGAL_SUPPORT),true)
CGAL_TESTS= \
	UnlabeledMatcher_test \
	MapAlignment_test
endif

ifeq ($(LIBSVM_SUPPORT),true)
LIBSVM_TESTS= \
	RTModel_test \
	RTPredict_test
endif

TOPP_TOOL_TESTS = ${UNCONDITIONAL_TESTS} ${LIBSVM_TESTS} ${CGAL_TESTS}

# Tests for processing pipelines using several tools

UNCONDITIONAL_PIPELINE_TESTS = \
	MyoQuant_test

ifeq ($(CGAL_SUPPORT),true)
CGAL_PIPELINE_TESTS = \
	PairwiseSuperposition_test
endif
TOPP_PIPELINE_TESTS = ${UNCONDITIONAL_PIPELINE_TESTS} ${CGAL_PIPELINE_TESTS}

# Files to remove
CLEANUP = ${RM} -f *.tmp *.log core.* *.TMP region* feature* model* plot.gp *.output

# Shortcut to TOPP bin path.
TOPP = ${OPENMS_PATH}/bin/

# Diff that allows for slight numeric differences
REL_DEVIATION = 1.05 # relative deviation allowed
ABS_DEVIATION = 0.01 # absolute deviation allowed
DEVIATION = $(REL_DEVIATION) $(ABS_DEVIATION)
NUMERICDIFF = NumericDiff

#######################################################################
# DEFAULT TARGET
# Here all TOPP tool tests are executed.
# Afterwards all temporary data files are removed (*.tmp)
default:  help $(TOPP_TOOL_TESTS) $(TOPP_PIPELINE_TESTS)
	@$(ECHO) ""
	@$(ECHO) "-------------------------------------------------------"
	@$(ECHO) ""
	@$(ECHO) "Congratulations!  All TOPP tests executed successfully."
	@$(ECHO) ""
	@ if test -z ${DEBUG}; then $(CLEANUP); fi

clean:
	@$(CLEANUP)

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
# DEBUG MODE
# Changes in debug mode are:
# - diff reports the differences (not only whether files are different)
# - temporary files (*.tmp) are not deleted after the execution
debug:
	@$(MAKE) DEBUG=1 default

# Diff that allows for slight numeric differences
DIFF_NUMERIC = \
	./$(NUMERICDIFF) $(1) $(2) $(3)

# Diff that ignores whitespaces and empty lines
DIFF_DEBUG = $(call DIFF_NUMERIC,$1,$2,$(DEVIATION))
DIFF_QUIET = $(call DIFF_NUMERIC,$1,$2,$(DEVIATION) -Q)

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
RUN_PROG_OPT_INTERNAL = echo "  Executed: $(TOPP)$(1) $(2)"; $(TOPP)$(1) $(2)
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

# 'diff' that allows for slight deviation of numbers
${NUMERICDIFF}:  NumericDiff.C

###########################################################################################
# Tests for individual TOPP tools.

AdditiveSeries_test: $(NUMERICDIFF)
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


AnalysisXMLMerger_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "AnalysisXMLMerger:"
	@$(call RUN_PROG_OPT,AnalysisXMLMerger,-ini AnalysisXMLMerger_1_parameters.ini)
	@$(call TEST_FILE_EQUAL,AnalysisXMLMerger_1_output.analysisXML,AnalysisXMLMerger_1_output.tmp)


BaselineFilter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "BaselineFilter:"
	@$(call RUN_PROG_OPT,BaselineFilter, -in RawData_input.mzData -out BaselineFilter_output.tmp -resampling 0.2 -struc_elem_length 1.5)
	@$(call TEST_FILE_EQUAL,BaselineFilter_output.tmp, BaselineFilter_output.mzData)


DBExporter_test: $(NUMERICDIFF)
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

DBImporter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "DBImporter:"
	@if test -f DBImporter_input.ini; then \
		$(ECHO) "- Test 1: -init"; \
		$(call RUN_PROG_OPT,DBImporter, -init -ini DBImporter_input.ini); \
		$(ECHO) "- Test 2: -in"; \
		$(call RUN_PROG_OPT,DBImporter, -in DBImporter_input.mzData -ini DBImporter_input.ini); \
		$(ECHO) "  The actual test is perfomed by DBExporter!"; \
	else \
		$(ECHO) "  Cannot run test as there is no INI file 'DBImporter_input.ini'"; \
	fi;


DTAExtractor_test: $(NUMERICDIFF)
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


MapDewarper_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "MapDewarper:"
	@$(call RUN_PROG_OPT,MapDewarper,-feat MapDewarper_1_feat.xml -grid MapDewarper_1_grid.xml -out MapDewarper_1_out.tmp )
	@$(call TEST_FILE_EQUAL,MapDewarper_1_out.xml,MapDewarper_1_out.tmp)

FeatureFinder_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "FeatureFinder:"

	@$(ECHO) "- Test 1:"
	@$(call RUN_PROG_OPT,FeatureFinder,\
	-ini FeatureFinder_1_parameters.ini \
	-in FeatureFinder_1_input.mzData \
	-out FeatureFinder_1_output.tmp)
	@$(call TEST_FILE_EQUAL,FeatureFinder_1_output.xml)


FeaturePairSplitter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "FeaturePairSplitter:"
	@$(call RUN_PROG_OPT,FeaturePairSplitter,\
	-in FeaturePairSplitter_1_pairs.xml \
	-out1 FeaturePairSplitter_1_out1.tmp \
	-out2 FeaturePairSplitter_1_out2.tmp \
	-qual FeaturePairSplitter_1_qual.tmp \
	-dump FeaturePairSplitter_1_dump.tmp)
	@$(MV) FeaturePairSplitter_1_dump.tmp.gp FeaturePairSplitter_1_dump_gp.tmp
	@$(SED) -i 's/^\#.*$$/# (all comments removed)/g' FeaturePairSplitter_1_dump_gp.tmp FeaturePairSplitter_1_dump.tmp
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_out1.xml)
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_out2.xml)
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_qual.wsv)
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_dump.wsv)
	@$(call TEST_FILE_EQUAL,FeaturePairSplitter_1_dump_gp.gp)


FileConverter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "FileConverter:"

	@$(ECHO) "- Test 1: MzData -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_1_input.txt -in_type MZDaTa -out FileConverter.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileConverter_1_output.dta2d,FileConverter.tmp)

	@$(ECHO) "- Test 2: DTA2D -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_2_input.txt -in_type DTA2d -out FileConverter.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileConverter_2_output.dta2d,FileConverter.tmp)

	@$(ECHO) "- Test 3: FeatureFile -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_3_input.txt -in_type FeatureFile -out FileConverter.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileConverter_3_output.dta2d,FileConverter.tmp)

	@$(ECHO) "- Test 4: MzXML -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_4_input.txt -in_type MZxml -out FileConverter.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileConverter_4_output.dta2d,FileConverter.tmp)

ifeq ($(ANDIMS_SUPPORT),true)
	@$(ECHO) "- Test 5: netCDF -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_5_input.txt -in_type ANDIMS -out FileConverter.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileConverter_5_output.dta2d,FileConverter.tmp)
endif
	@$(ECHO) "- Test 6: DTA2D -> MzXML -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_6_input.dta2d -out FileConverter.tmp -out_type mzxml)
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter.tmp -in_type mzxml -out FileConverter2.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileConverter_6_input.dta2d,FileConverter2.tmp)

	@$(ECHO) "- Test 7: DTA2D -> MzData -> DTA2D"
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter_6_input.dta2d -out FileConverter.tmp -out_type mzData)
	@$(call RUN_PROG_OPT,FileConverter,-in FileConverter.tmp -in_type mzdata -out FileConverter2.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileConverter_6_input.dta2d,FileConverter2.tmp)


FileFilter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "FileFilter:"

	@$(ECHO) "- Test 1: RT, m/z test, int (upper bound)"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter_1_out.tmp -rt :30 -mz :1000 -int :20000)
	@$(call TEST_FILE_EQUAL,FileFilter_1_output.mzData,FileFilter_1_out.tmp)

	@$(ECHO) "- Test 2: RT, m/z test, int (lower bound)"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter_2_out.tmp -rt 30: -mz 1000: -int 100:)
	@$(call TEST_FILE_EQUAL,FileFilter_2_output.mzData,FileFilter_2_out.tmp)

	@$(ECHO) "- Test 3: level"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_1_input.mzData -out FileFilter_3_out.tmp -level 2)
	@$(call TEST_FILE_EQUAL,FileFilter_3_output.mzData,FileFilter_3_out.tmp)

	@$(ECHO) "- Test 4: remove_zoom"
	@$(call RUN_PROG_OPT,FileFilter,-in FileFilter_4_input.mzData -out FileFilter_4_out.tmp -remove_zoom)
	@$(call TEST_FILE_EQUAL,FileFilter_4_output.mzData,FileFilter_4_out.tmp)

FileInfo_test: $(NUMERICDIFF)
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

ifeq ($(ANDIMS_SUPPORT),true)
	@$(ECHO) "- Test 6: netCDF file test (with meta info)"
	@$(call RUN_PROG_OPT,FileInfo,-in FileInfo_6_input.txt -in_type ANDIMS -m > FileInfo.tmp)
	@$(call TEST_FILE_EQUAL,FileInfo_6_output.txt,FileInfo.tmp)
endif

FileMerger_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "FileMerger:"

	@$(ECHO) "- Test 1: 3 Files (-rt_file)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_1_ini.txt -out FileMerger.tmp -rt_file)
	@$(call RUN_PROG_OPT,FileConverter,-in FileMerger.tmp -in_type MzData -out FileMerger2.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileMerger_1_output.txt,FileMerger2.tmp)

	@$(ECHO) "- Test 2: 2 Files (-rt_file -in_type)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_2_ini.txt -in_type DTA -out FileMerger.tmp -rt_file)
	@$(call RUN_PROG_OPT,FileConverter,-in FileMerger.tmp -in_type MzData -out FileMerger2.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileMerger_2_output.txt,FileMerger2.tmp)

	@$(ECHO) "- Test 3: 2 Files (-rt_auto)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_3_ini.txt -out FileMerger.tmp -rt_auto)
	@$(call RUN_PROG_OPT,FileConverter,-in FileMerger.tmp -in_type MzData -out FileMerger2.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileMerger_3_output.txt,FileMerger2.tmp)

	@$(ECHO) "- Test 4: 2 Files (several files)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_4_ini.txt -out FileMerger.tmp)
	@$(call RUN_PROG_OPT,FileConverter,-in FileMerger.tmp -in_type MzData -out FileMerger2.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,FileMerger_4_output.txt,FileMerger2.tmp)

	@$(ECHO) "- Test 5: 3 Files (guess rt from filename, -rt_from_filename, -user_ms_level, -ms_level)"
	@$(call RUN_PROG_OPT,FileMerger,-file_list FileMerger_5_ini.txt -rt_from_filename -user_ms_level -ms_level 2  -out FileMerger5.tmp)
	@$(call TEST_FILE_EQUAL,FileMerger_5_output.txt,FileMerger5.tmp)

IDFilter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "IDFilter:"

	@$(ECHO) "- Test 1: sequence file"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_1_input.analysisXML -out IDFilter_1_output.tmp -sequences_file IDFilter_1_input.fas)
	@$(call TEST_FILE_EQUAL,IDFilter_1_output.analysisXML,IDFilter_1_output.tmp)

	@$(ECHO) "- Test 2: threshold fractions"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_1_input.analysisXML -out IDFilter_2_output.tmp -pep_fraction 1 -prot_fraction 1)
	@$(call TEST_FILE_EQUAL,IDFilter_2_output.analysisXML,IDFilter_2_output.tmp)

	@$(ECHO) "- Test 3: exclusion peptides"
	@$(call RUN_PROG_OPT,IDFilter,-in IDFilter_1_input.analysisXML -out IDFilter_3_output.tmp -exclusion_peptides_file IDFilter_3_input.analysisXML)
	@$(call TEST_FILE_EQUAL,IDFilter_3_output.analysisXML,IDFilter_3_output.tmp)

LabeledMatcher_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "LabeledMatcher:"
	@$(call RUN_PROG_OPT,LabeledMatcher,-in LabeledMatcher_1_input.feat -ini LabeledMatcher_1_parameters.ini -out LabeledMatcher_1_output.tmp)
	@$(call TEST_FILE_EQUAL,LabeledMatcher_1_output.pairs,LabeledMatcher_1_output.tmp)


MapAlignment_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "MapAlignment:"
	
	@$(ECHO) "- Test 1: feature files"
	@$(call RUN_PROG_OPT,MapAlignment, -ini MapAlignment_1_feature_parameters.ini -out MapAlignment_1_feature_output.tmp)
	@$(call TEST_FILE_EQUAL,MapAlignment_1_feature_consensus.xml, MapAlignment_1_feature_output.tmp)

	@$(ECHO) "- Test 2: peak files"
	@$(call RUN_PROG_OPT,MapAlignment, -ini MapAlignment_1_peak_parameters.ini -out MapAlignment_1_peak_output.tmp)
	@$(call TEST_FILE_EQUAL,MapAlignment_1_peak_consensus.xml, MapAlignment_1_peak_output.tmp)

	@$(ECHO) "- Test 3: consensus files"
	@$(call RUN_PROG_OPT,MapAlignment, -ini MapAlignment_1_consensus_parameters.ini -out MapAlignment_1_consensus_output.tmp)
	@$(call TEST_FILE_EQUAL,MapAlignment_1_consensus_consensus.xml, MapAlignment_1_consensus_output.tmp)


MapMatcher_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "MapMatcher:"
	@$(call RUN_PROG_OPT,MapMatcher,-pairs MapMatcher_1_input.pairs -grid MapMatcher_1_input.grid -out MapMatcher_1_output_grid.tmp)
	@$(call TEST_FILE_EQUAL,MapMatcher_1_output.grid,MapMatcher_1_output_grid.tmp)

MapStatistics_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "MapStatistics:"
	@$(call RUN_PROG_OPT,MapStatistics,-in MapStatistics_rawdata.mzData -out MapStatistics_rawdata_results.tmp)
	@$(call TEST_FILE_EQUAL,MapStatistics_rawdata_results.tmp, MapStatistics_rawdata_results.xml)
	@$(call RUN_PROG_OPT,MapStatistics,-in MapStatistics_featdata.feat -out MapStatistics_featdata_results.tmp)
	@$(call TEST_FILE_EQUAL,MapStatistics_featdata_results.tmp, MapStatistics_featdata_results.xml)

MascotAdapter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "MascotAdapter:"
	@$(call RUN_PROG_OPT,MascotAdapter,-ini MascotAdapter_1_parameters.ini -mascot_in)
	@$(call TEST_FILE_EQUAL,MascotAdapter_1_output.mascot_in,MascotAdapter_1_output.tmp)
	@$(call RUN_PROG_OPT,MascotAdapter,-in MascotAdapter_2_input.mascotXML -out MascotAdapter_2_output.tmp -mascot_out)
	@$(call TEST_FILE_EQUAL,MascotAdapter_2_output.analysisXML,MascotAdapter_2_output.tmp)


NoiseFilter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "NoiseFilter:"

	@$(ECHO) "- Test 1: Gaussian filter"
	@$(call RUN_PROG_OPT,NoiseFilter, -ini NoiseFilter_1_parameters.ini -filter_type gaussian -in  RawData_input.mzData -out NoiseFilter_1_output.tmp)
	@$(call TEST_FILE_EQUAL,NoiseFilter_1_output.tmp, NoiseFilter_1_output.mzData)

	@$(ECHO) "- Test 2: Savitzky Golay filter"
	@$(call RUN_PROG_OPT,NoiseFilter,	-ini NoiseFilter_2_parameters.ini -filter_type sgolay -resampling 0.2 -in  RawData_input.mzData -out NoiseFilter_2_output.tmp)
	@$(call TEST_FILE_EQUAL,NoiseFilter_2_output.tmp, NoiseFilter_2_output.mzData)


PeakPicker_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "PeakPicker:"
	
	@$(ECHO) "-Test1 :PeakPicker without deconvolution" 	 
	@$(call RUN_PROG_OPT,PeakPicker, -ini PeakPicker_parameters.ini -in RawData_input.mzData -out PeakPicker_output.tmp)
	@$(call TEST_FILE_EQUAL,PeakPicker_output.tmp, PeakPicker_output.mzData)

	@$(ECHO) ""
	@$(ECHO) "-Test2 :PeakPicker with deconvolution"
	@$(call RUN_PROG_OPT,PeakPicker, -ini PeakPicker_deconv_parameters.ini -in PeakPicker_deconv_input.mzData -out PeakPicker_deconv_output.tmp)
	@$(call TEST_FILE_EQUAL,PeakPicker_deconv_output.tmp, PeakPicker_deconv_output.mzData)


Resampler_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "Resampler:"

	@$(ECHO) "- Test 1: output mzData"
	@$(call RUN_PROG_OPT,Resampler,-in Resampler_1_input.mzData -out Resampler_1_output.tmp -mz 685:697 -cols_mz 20 -rt 1230:1260 -rows_rt 15)
	@$(call RUN_PROG_OPT,FileConverter,-in Resampler_1_output.tmp -in_type mzData -out Resampler_1_output_dta2d.tmp -out_type dta2d)
	@$(call TEST_FILE_EQUAL,Resampler_1_output.dta2d,Resampler_1_output_dta2d.tmp)

	@$(ECHO) "- Test 2: output pgm (using default settings)"
	@$(call RUN_PROG_OPT,Resampler,-in Resampler_1_input.mzData -pgm Resampler_2_output.tmp -cols_mz 20 -rows_rt 15)
	@$(SED) -i 's/^\# generated by TOPP Resampler on .*$$/# generated by TOPP Resampler on <date> <time>/g' Resampler_2_output.tmp
	@$(call TEST_FILE_EQUAL,Resampler_2_output.pgm)

	@$(ECHO) "- Test 3: output pgm (using all bells & whistles)"
	@$(call RUN_PROG_OPT,Resampler,-in Resampler_1_input.mzData  -pgm Resampler_3_output.tmp -rt : -mz : -cols_mz 15 -rows_rt 20 -transpose -maxval 230 -reverse -gamma .2 -scale 0.0015)
	@$(SED) -i 's/^\# generated by TOPP Resampler on .*$$/# generated by TOPP Resampler on <date> <time>/g' Resampler_3_output.tmp
	@$(call TEST_FILE_EQUAL,Resampler_3_output.pgm)

RTModel_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "RTModel:"
	@$(call RUN_PROG_OPT,RTModel,-in RTModel_1_input.analysisXML -out RTModel_1_output.tmp -ini RTModel_1_parameters.ini)
	@$(SED) -i 's/^probA.*$$/# (removed probA line)/g' RTModel_1_output.tmp RTModel_1_output.tmp
	@$(call TEST_FILE_EQUAL,RTModel_1_output.tmp,RTModel_1_output.model)


RTPredict_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "RTPredict:"
	@$(call RUN_PROG_OPT,RTPredict,-in RTPredict_1_input.analysisXML -out RTPredict_1_output.tmp -total_gradient_time 3000 -svm_model RTPredict_1_input.model)
	@$(call TEST_FILE_EQUAL,RTPredict_1_output.tmp,RTPredict_1_output.analysisXML)

InspectAdapter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "InspectAdapter:"
	
	@$(ECHO) "- Test 1: -inspect_in"
	@$(call RUN_PROG_OPT,InspectAdapter,-ini InspectAdapter_1_parameters.ini -inspect_in)
	@$(SED) -i 's/^spectra,.*$$/spectra,dummy.mzXML/g' InspectAdapter_1_output.tmp
	@$(SED) -i 's/^db,.*$$/db,dummy.trie/g' InspectAdapter_1_output.tmp
	@$(call TEST_FILE_EQUAL,InspectAdapter_1_output.inspect_in,InspectAdapter_1_output.tmp)
	
	@$(ECHO) "- Test 2: -inspect_out"
	@$(call RUN_PROG_OPT,InspectAdapter,-ini InspectAdapter_1_parameters.ini -in InspectAdapter.out -out InspectAdapter_2_output.tmp -inspect_out)
	@$(SED) -i 's/^.*date_group_0.*$$/# date group removed/g' InspectAdapter_2_output.tmp
	@$(call TEST_FILE_EQUAL,InspectAdapter_2_output.analysisXML,InspectAdapter_2_output.tmp)
	@$(RM) -rf Inspect_FASTAFile_test4.trie Inspect_FASTAFile_test4.index

SequestAdapter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "SequestAdapter:"
	
	@$(RM) -rf ./tmp/*.dta_*  ./tmp/sequest_run.bat ./tmp/temp.sequest.in
	@$(ECHO) "- Test 1: -sequest_in"
	@$(call RUN_PROG_OPT,SequestAdapter,-ini SequestAdapter_1_parameters.ini -in Inspect_Sequest.mzXML -out SequestAdapter_1_output.tmp -sequest_in)
	@$(call TEST_FILE_EQUAL,SequestAdapter_1_output.sequest_in,SequestAdapter_1_output.tmp)
	@$(RM) -rf Inspect_Sequest*.dta_*
	
	@$(ECHO) "- Test 2: -sequest_out"
	@$(call RUN_PROG_OPT,SequestAdapter,-ini SequestAdapter_2_parameters.ini -in tmp -out SequestAdapter_2_output.tmp -sequest_out)
	@$(SED) -i 's/^.*date_group_1.*$$/# date group removed/g' SequestAdapter_2_output.tmp
	@$(call TEST_FILE_EQUAL,SequestAdapter_2_output.analysisXML,SequestAdapter_2_output.tmp)
	@$(RM) -rf ./tmp/sequest_run.bat ./tmp/temp.sequest.in

SpectraFilter_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "SpectraFilter:"
	
	@$(ECHO) "- Test 1: -filter SqrtMower"
	@$(call RUN_PROG_OPT,SpectraFilter,-in SpectraFilter_1_input.mzData -out SpectraFilter_1_output.mzData.tmp -filters SqrtMower)
	@$(call TEST_FILE_EQUAL,SpectraFilter_1_output.mzData,SpectraFilter_1_output.mzData.tmp)
	
	@$(ECHO) "- Test 2: -filter WindowMower"
	@$(call RUN_PROG_OPT,SpectraFilter,-in SpectraFilter_2_input.mzData -out SpectraFilter_2_output.mzData.tmp -filters WindowMower)
	@$(call TEST_FILE_EQUAL,SpectraFilter_2_output.mzData,SpectraFilter_2_output.mzData.tmp)
	
	@$(ECHO) "- Test 3: -filter WindowMower -ini SpectraFilter_3_parameters.ini"
	@$(call RUN_PROG_OPT,SpectraFilter,-in SpectraFilter_3_input.mzData -out SpectraFilter_3_output.mzData.tmp -ini  SpectraFilter_3_parameters.ini -filters WindowMower)
	@$(call TEST_FILE_EQUAL,SpectraFilter_3_output.mzData,SpectraFilter_3_output.mzData.tmp)        

UnlabeledMatcher_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "UnlabeledMatcher:"
	@$(call RUN_PROG_OPT,UnlabeledMatcher,-ini UnlabeledMatcher_1_parameters.ini)
	@$(call TEST_FILE_EQUAL,UnlabeledMatcher_1_pairs.xml)
	@$(call TEST_FILE_EQUAL,UnlabeledMatcher_1_grid.xml)


###########################################################################################
## Tests for processing pipelines.

MyoQuant_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "MyoQuant_test:"

	@$(ECHO) "- Test 1: Additive Series using intensity ratio of tryptic peptides T11hu / T10ho, no dewarping"

	@$(foreach number, 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 , \
	  ${ECHO} Finding features in LC-MS map ${number}; \
	  $(call RUN_PROG_OPT,FeatureFinder,-ini MyoQuant_1_parameters.ini -in MyoQuant_1_input${number}.mzData -out MyoQuant_1_feat${number}_xml.tmp); \
	  )
	@$(ECHO) "	AdditiveSeries, linear regression";
	@$(call RUN_PROG_OPT,AdditiveSeries,-ini MyoQuant_1_parameters.ini)
	@$(call TEST_FILE_EQUAL,MyoQuant_1_output.xml)

	@$(ECHO) ""
	@$(ECHO) "- Test 2: Additive Series using intensity ratio of tryptic peptides T11hu / T10ho, with dewarping"
	@$(ECHO) "	... todo ..."

PairwiseSuperposition_test: $(NUMERICDIFF)
	@$(ECHO) ""
	@$(ECHO) "PairwiseSuperposition_test:"
	@$(ECHO) "- Test 1: UnlabeledMatcher, MapMatcher, and MapDewarper"
	@$(ECHO) "	PairwiseSuperposition_1_input1.feat shall be transformed into PairwiseSuperposition_1_output1.feat"
	@$(ECHO) "	which is much closer to PairwiseSuperposition_1_input2.feat"
	@$(ECHO) "  Find pairs using simple UnlabeledMatcher..."
	@$(call RUN_PROG_OPT,UnlabeledMatcher,-ini PairwiseSuperposition_1_parameters.ini -in1 PairwiseSuperposition_1_input1.feat -in2 PairwiseSuperposition_1_input2.feat -pairs PairwiseSuperposition_1_output1_pairs.tmp -grid PairwiseSuperposition_1_output1_grid.tmp)
	@$(call TEST_FILE_EQUAL,PairwiseSuperposition_1_output1.grid,PairwiseSuperposition_1_output1_grid.tmp)
	@$(call TEST_FILE_EQUAL,PairwiseSuperposition_1_output1.pairs,PairwiseSuperposition_1_output1_pairs.tmp)
	@$(ECHO) "  Find affine map by linear regression using MapMatcher..."
	@$(call RUN_PROG_OPT,MapMatcher,-pairs PairwiseSuperposition_1_output1.pairs -grid PairwiseSuperposition_1_output1.grid -out PairwiseSuperposition_1_output2_grid.tmp)
	@$(call TEST_FILE_EQUAL,PairwiseSuperposition_1_output2.grid,PairwiseSuperposition_1_output2_grid.tmp)
	@$(ECHO) "  Apply affine map using MapDewarper..."
	@$(call RUN_PROG_OPT,MapDewarper,-feat PairwiseSuperposition_1_input1.feat -grid PairwiseSuperposition_1_output2.grid -out PairwiseSuperposition_1_output1_feat.tmp)
	@$(call TEST_FILE_EQUAL,PairwiseSuperposition_1_output1.feat,PairwiseSuperposition_1_output1_feat.tmp)
