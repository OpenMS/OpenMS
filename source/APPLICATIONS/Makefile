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
# $Id: Makefile,v 1.13 2006/06/07 06:13:36 marc_sturm Exp $
# $Author: marc_sturm $
# $Maintainer:  $
# --------------------------------------------------------------------------
#

include ../config.mak
-include .Dependencies

DIRECTORY=APPLICATIONS

THISLIB=$(LIBNAME)

CPP_SOURCES= \

DIRS=\
	CLUSTERANALYSIS\
	TOPP

#FEATUREFINDER
#no makefile in FEATUREFINDER dir => make clean fails

default:	dirs

dirs:	dummy
	@$(MAKE) subdirs SUBDIRS="$(DIRS)" TARGET="all"

depend:	dummy
	@$(MAKE) subdirs SUBDIRS="$(DIRS)" TARGET=depend

subdirs:	dummy
	@set -e; for i in $(SUBDIRS); do cd $$i; if $(MAKE) $(TARGET); then cd .. ; else cd .. ; exit 1 ; fi ; done

default_clean:
	@-$(RM) -rf .Dependencies *~  core.* $(OBJECTS) $(MOC_SOURCES) $(ADD_MOC_SOURCES) $(EXECUTABLES) $(EXECUTABLE)
	@$(MAKE) subdirs SUBDIRS="$(DIRS)" TARGET=clean

collectlib:

dummy:




