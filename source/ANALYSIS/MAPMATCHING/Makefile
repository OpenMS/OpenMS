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
# $Maintainer: Eva Lange $
# --------------------------------------------------------------------------

DIRECTORY=ANALYSIS/MAPMATCHING

THISLIB=$(LIBNAME)

include ../../config_defs.mak

ifeq ($(CGAL_SUPPORT),true)
CGAL_SOURCES=\
	BaseAlignment.C\
	StarAlignment.C\
	DelaunayPairFinder.C
endif

CPP_SOURCES=\
  DBaseMapMatcher.C\
	DBaseMapping.C\
	DFeaturePair.C\
	DFeaturePairVector.C\
	DGridCell.C\
	DGrid.C\
	DMapDewarper.C\
	DLinearMapping.C\
	PairMatcher.C\
	BasePairwiseMapMatcher.C\
	PoseClusteringPairwiseMapMatcher.C\
	BaseSuperimposer.C\
	PoseClusteringShiftSuperimposer.C\
	PoseClusteringAffineSuperimposer.C\
	BasePairFinder.C\
	SimplePairFinder.C\
	Group.C\
	IndexTuple.C\
	DMapMatcherRegression.C

ADD_CPP_SOURCES = $(CGAL_SOURCES)

include ../../common.mak
-include .Dependencies
