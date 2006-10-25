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

# tool tests in alphabetical order
TOPP_TESTS= \
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


# DEFAULT TARGET
# Here all TOPP tool tests are executed. Afterwards all temporary data files are removed (*.tmp)
default:  $(TOPP_TESTS)
	@${ECHO} ""
	@${ECHO} "All TOPP tests executed successfully"
	@${ECHO} ""
	@ if test -z ${DEBUG}; then ${RM} -f *.tmp; fi

# DEBUG MODE
# You can activate the debug mode with 'make debug'.
# Changes in debug mode are:
# - diff reports the differences (not only whether files are different)
# - temporary files (*.tmp) are not deleted after the execution
debug:
	@${MAKE} DEBUG=1 default


# Diff that ignores whitespaces and empty lines
ifeq ("$(DEBUG)","")
  DIFF_WITH_OPT = ${DIFF} -bB --brief
else
  DIFF_WITH_OPT = ${DIFF} -bBu
endif

# Shortcut to TOPP bin path
TOPP = ${OPENMS_PATH}/bin/


##  Tool tests coming up soon...
#
##  Marc (?) :
# DBImporter_test
# DBExporter_test
# LabeledMatcher_test
#
##  Clemens:
# FeatureFinder_test
# UnlabeledMatcher_test
# AdditiveSeries_test


AdditiveSeries_test:
	@${ECHO} ""
	@${ECHO} "AdditiveSeries:"


AnalysisXMLMerger_test:
	@${ECHO} ""
	@${ECHO} "AnalysisXMLMerger:"


BaselineFilter_test:
	@${ECHO} ""
	@${ECHO} "BaselineFilter:"


DBExporter_test:
	@${ECHO} ""
	@${ECHO} "DBExporter:"


DBImporter_test:
	@${ECHO} ""
	@${ECHO} "DBImporter:"


DTAExtractor_test:
	@${ECHO} ""
	@${ECHO} "DTAExtractor:"

	@${ECHO} "- Test 1: rt"
	@${TOPP}DTAExtractor -in DTAExtractor_1_input.mzData -out DTAExtractor -rt :61
	@${MV} DTAExtractor_RT60.dta DTAExtractor.tmp
	@${DIFF_WITH_OPT} DTAExtractor.tmp DTAExtractor_1_output.dta

	@${ECHO} "- Test 2: level"
	@${TOPP}DTAExtractor -in DTAExtractor_1_input.mzData -out DTAExtractor -level 1
	@${MV} DTAExtractor_RT60.dta DTAExtractor.tmp
	@${DIFF_WITH_OPT} DTAExtractor.tmp DTAExtractor_2_output.dta

	@${ECHO} "- Test 3: precursor mz, level"
	@${TOPP}DTAExtractor -in DTAExtractor_1_input.mzData -out DTAExtractor -level 2 -mz :1000
	@${MV} DTAExtractor_RT140_MZ5.dta DTAExtractor.tmp
	@${DIFF_WITH_OPT} DTAExtractor.tmp DTAExtractor_3_output.dta


Dewarper_test:
	@${ECHO} ""
	@${ECHO} "Dewarper:"


FeatureFinder_test:
	@${ECHO} ""
	@${ECHO} "FeatureFinder:"


FeaturePairSplitter_test:
	@${ECHO} ""
	@${ECHO} "FeaturePairSplitter:"


FileConverter_test:
	@${ECHO} ""
	@${ECHO} "FileConverter:"


FileFilter_test:
	@${ECHO} ""
	@${ECHO} "FileFilter:"

	@${ECHO} "- Test 1: RT, m/z test, int (upper bound)"
	@${TOPP}FileFilter -in FileFilter_1_input.mzData -out FileFilter.tmp -rt :30 -mz :1000 -int :20000
	@${DIFF_WITH_OPT} FileFilter.tmp FileFilter_1_output.mzData

	@${ECHO} "- Test 2: RT, m/z test, int (lower bound)"
	@${TOPP}FileFilter -in FileFilter_1_input.mzData -out FileFilter.tmp -rt 30: -mz 1000: -int 100:
	@${DIFF_WITH_OPT} FileFilter.tmp FileFilter_2_output.mzData

	@${ECHO} "- Test 3: level"
	@${TOPP}FileFilter -in FileFilter_1_input.mzData -out FileFilter.tmp -level 2
	@${DIFF_WITH_OPT} FileFilter.tmp FileFilter_3_output.mzData

	@${ECHO} "- Test 4: remove_zoom"
	@${TOPP}FileFilter -in FileFilter_4_input.mzData -out FileFilter.tmp -remove_zoom
	@${DIFF_WITH_OPT} FileFilter.tmp FileFilter_4_output.mzData


FileInfo_test:
	@${ECHO} ""
	@${ECHO} "FileInfo:"

	@${ECHO} "- Test 1: DTA file test"
	@${TOPP}FileInfo -in FileInfo_1_input.txt -in_type DTa > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_1_output.txt

	@${ECHO} "- Test 2: DTA2D file test"
	@${TOPP}FileInfo -in FileInfo_2_input.dta2d > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_2_output.txt

	@${ECHO} "- Test 3: feature file test (with meta info)"
	@${TOPP}FileInfo -in FileInfo_3_input.feat -m > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_3_output.txt

	@${ECHO} "- Test 4: MzXML file test (with meta info)"
	@${TOPP}FileInfo -in FileInfo_4_input.mzXML -m > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_4_output.txt

	@${ECHO} "- Test 5: MzData file test (with meta info)"
	@${TOPP}FileInfo -in FileInfo_5_input.txt -in_type MzDaTa -m > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_5_output.txt


FileMerger_test:
	@${ECHO} ""
	@${ECHO} "FileMerger:"


IDFilter_test:
	@${ECHO} ""
	@${ECHO} "IDFilter:"


InspectAdapter_test:
	@${ECHO} ""
	@${ECHO} "InspectAdapter:"


LabeledMatcher_test:
	@${ECHO} ""
	@${ECHO} "LabeledMatcher:"


MapAlignment_test:
	@${ECHO} ""
	@${ECHO} "MapAlignment:"


MapMatcher_test:
	@${ECHO} ""
	@${ECHO} "MapMatcher:"


MapStatistics_test:
	@${ECHO} ""
	@${ECHO} "MapStatistics:"


MascotAdapter_test:
	@${ECHO} ""
	@${ECHO} "MascotAdapter:"


NoiseFilter_test:
	@${ECHO} ""
	@${ECHO} "NoiseFilter:"


PeakPicker_test:
	@${ECHO} ""
	@${ECHO} "PeakPicker:"


RTModel_test:
	@${ECHO} ""
	@${ECHO} "RTModel:"


RTPredict_test:
	@${ECHO} ""
	@${ECHO} "RTPredict:"


SequestAdapter_test:
	@${ECHO} ""
	@${ECHO} "SequestAdapter:"


SpectraFilter_test:
	@${ECHO} ""
	@${ECHO} "SpectraFilter:"


TOPPView_test:
	@${ECHO} ""
	@${ECHO} "TOPPView:"


UnlabeledMatcher_test:
	@${ECHO} ""
	@${ECHO} "UnlabeledMatcher:"


