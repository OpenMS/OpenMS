# -*- Mode: make; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2004 -- Oliver Kohlbacher, Knut Reinert
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
# $Id: Makefile,v 1.5 2006/03/10 14:55:27 elange Exp $
# $Author: elange $
# $Maintainer:  $
# --------------------------------------------------------------------------
#

DIRECTORY=FILTERING/BASELINE

THISLIB=$(LIBNAME)

include ../../config_defs.mak

CPP_SOURCES= \
#conditional compiles
ifeq ($(GSL_SUPPORT),true)
	ADD_CPP_SOURCES=\
	DMorphFilter.C\
	DTopHatFilter.C
endif

include ../../common.mak
-include .Dependencies

filter: filter.C 
	${CXX} -g $(LDFLAGS) $(OPENMS_INCLUDES) filter.C $(LIBS) -L/home/arabidopsis/software/Linux/OpenMS/qt/lib -lqt-mt -o filter

filter_2D: filter_2D.C 
	${CXX} -g $(LDFLAGS) $(OPENMS_INCLUDES) filter_2D.C $(LIBS) -L/home/arabidopsis2/software/Linux/OpenMS/qt/lib -lqt-mt -o filter_2D

LIBS=$(GUI_LIBS)  

filt: filt.C 
	${CXX} -g $(LDFLAGS) $(OPENMS_INCLUDES) filt.C $(LIBS) -L/home/arabidopsis/software/Linux/OpenMS/qt/lib -lqt-mt -o filt

TopHatFiltering: TopHatFiltering.C 
	${CXX} -g $(LDFLAGS) $(OPENMS_INCLUDES) TopHatFiltering.C $(LIBS) -L/home/arabidopsis/software/Linux/OpenMS/qt/lib -lqt-mt -o TopHatFiltering
