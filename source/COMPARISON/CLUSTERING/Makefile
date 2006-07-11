# -*- Mode: make; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Ulf Leser, Knut Reinert
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


DIRECTORY=COMPARISON/CLUSTERING

THISLIB=$(LIBNAME)

include ../../config_defs.mak


#conditional compiles
ifeq ($(ENABLE_GUI),true)
  DB_SOURCES=\
  ClusterIdentificator.C\
  Cluster.C\
  BinnedRep.C \
  ClusterFactory.C\
  ClusterExperiment.C\
  ClusterRun.C\
  Analysis.C\
  ClusterRunAnalysisLess.C\
	AnalysisFunctor.C\
  ClusterSpectrum.C\
	ClusterCompareFunctor.C\
  ClusterFunctor.C\
  LinkageCluster.C\
  DistanceAnalyzer.C

  
  CPP_SOURCES=\
    ClusterNode.C\
    helper.C\
    SpectrumGenerator.C\
		ClusterExperimentXMLHandler.C\
    $(DB_SOURCES)
 
endif

  
include ../../common.mak
-include .Dependencies
OPENMS_INCLUDES=$(GUI_INCLUDES)
