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
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include ../../config_defs.mak

DIRECTORY=TEST/TOPP

#diff that ignores whitespaces and empty lines
DIFF_WITH_OPT = ${DIFF} -bB --brief

# default target 
# Here all TOPP tool tests are executed. Afterwards all temporary data files are removed (*.tmp)
# You can disable removing the tmp files by calling 'make DEGUB=1 <target>'
default: FileInfo_test FileFilter_test
	@${ECHO} 
	@${ECHO} All TOPP tests executed successfully 
	@${ECHO} 
	@ if test -z ${DEBUG}; then ${RM} -f *.tmp; fi 

FileInfo_test:
	@${ECHO}
	@${ECHO} FileInfo:
	
	@${ECHO} - Test 1: DTA file test
	@${OPENMS_PATH}/bin/FileInfo -in FileInfo_1_input.txt -in_type DTa > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_1_output.txt

	@${ECHO} - Test 2: DTA2D file test
	@${OPENMS_PATH}/bin/FileInfo -in FileInfo_2_input.dta2d > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_2_output.txt

	@${ECHO} - Test 3: feature file test
	@${OPENMS_PATH}/bin/FileInfo -in FileInfo_3_input.feat -m > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_3_output.txt

	@${ECHO} - Test 4: MzXML file test
	@${OPENMS_PATH}/bin/FileInfo -in FileInfo_4_input.mzXML -m > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_4_output.txt

	@${ECHO} - Test 5: MzData file test
	@${OPENMS_PATH}/bin/FileInfo -in FileInfo_5_input.txt -in_type MzDaTa -m > FileInfo.tmp
	@${DIFF_WITH_OPT} FileInfo.tmp FileInfo_5_output.txt


FileFilter_test:
	@${ECHO}
	@${ECHO} FileFilter:

	@${ECHO} - Test 1: RT, m/z test, int - upper bound
	@${OPENMS_PATH}/bin/FileFilter -in FileFilter_1_input.mzData -out FileFilter.tmp -rt :30 -mz :1000 -int :20000
	@${DIFF_WITH_OPT} FileFilter.tmp FileFilter_1_output.mzData

	@${ECHO} - Test 2: RT, m/z test, int - lower bound
	@${OPENMS_PATH}/bin/FileFilter -in FileFilter_1_input.mzData -out FileFilter.tmp -rt 30: -mz 1000: -int 100:
	@${DIFF_WITH_OPT} FileFilter.tmp FileFilter_2_output.mzData

	@${ECHO} - Test 3: level
	@${OPENMS_PATH}/bin/FileFilter -in FileFilter_1_input.mzData -out FileFilter.tmp -level 2
	@${DIFF_WITH_OPT} FileFilter.tmp FileFilter_3_output.mzData

	@${ECHO} - Test 4: remove_zoom
	@${OPENMS_PATH}/bin/FileFilter -in FileFilter_4_input.mzData -out FileFilter.tmp -remove_zoom
	@${DIFF_WITH_OPT} FileFilter.tmp FileFilter_4_output.mzData


#DTAExtractor
#DBImporter
#FileConverter
#DBExporter
#FileMerger
#LabeledMatcher