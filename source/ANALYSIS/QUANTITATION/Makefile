# -*- mode: makefile; tab-width: 2; -*-
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

DIRECTORY=ANALYSIS/QUANTITATION


include ../../config_defs.mak

CPP_SOURCES=\
		ItraqConstants.C\
  	ItraqChannelExtractor.C\
		ItraqQuantifier.C\
		ProteinInference.C


include ../../common.mak
-include .Dependencies

# For debugging the ID classes, create the file 'Makefile_debug' and 
# set 'ADD_CXXFLAGS' one or more of the following flags
#  -DDEBUG_ID_MAPPING  
#  -DDEBUG_ID_CONSENSUS
-include Makefile_debug
