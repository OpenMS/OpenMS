# -*- Mode: make; tab-width: 2; -*-
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
# $Maintainer:  $
# --------------------------------------------------------------------------
#

DIRECTORY=FORMAT/HANDLERS

THISLIB=$(LIBNAME)

include ../../config_defs.mak

CPP_SOURCES=\
	XMLHandler.C\
	SchemaHandler.C\
	MzDataHandler.C\
	MzXMLHandler.C\
	ANDIHandler.C\
	ParamXMLHandler.C\
	DFeatureMapHandler.C\
	AnalysisXMLHandler.C\
	DFeaturePairsHandler.C\
	DGridHandler.C\
	MzDataExpSettHandler.C\
	MascotXMLHandler.C\
	XMLSchemes.C\
	OMSSAXMLHandler.C\
	PTMXMLHandler.C\
	IsotopeXMLHandler.C\

ifeq ($(ANDIMS_SUPPORT),true)
	ifeq ($(NETCDF_SUPPORT),true)
		ANDI_SOURCES=\
			ANDIHandler.C
	endif
endif

ifeq ($(CGAL_SUPPORT),true)
	CGAL_SOURCES=\
		ConsensusXMLHandler.C
endif

ADD_CPP_SOURCES = $(ANDI_SOURCES) $(CGAL_SOURCES)

include ../../common.mak
-include .Dependencies





