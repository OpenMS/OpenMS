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
# $Maintainer:  $
# --------------------------------------------------------------------------
#

DIRECTORY=FORMAT/HANDLERS


include ../../config_defs.mak

CPP_SOURCES=\
	XMLHandler.C\
	MzDataHandler.C\
	MzXMLHandler.C\
	MzMLHandler.C\
	ParamXMLHandler.C\
	FeatureXMLHandler.C\
	MzDataExpSettHandler.C\
	MascotXMLHandler.C\
	PTMXMLHandler.C\
	ConsensusXMLHandler.C\
	XTandemInfileXMLHandler.C\
	UnimodXMLHandler.C

ifeq ($(ANDIMS_SUPPORT),true)
	ANDI_SOURCES=\
		ANDIHandler.C
endif

ADD_CPP_SOURCES = $(ANDI_SOURCES)

include ../../common.mak
-include .Dependencies





