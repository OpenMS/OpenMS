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
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

all:	default

include ../config_defs.mak

DIRECTORY=ANALYSIS

THISLIB=$(LIBNAME)

-include .Dependencies

# Take care of subdirectories

SECOND_LEVEL_DIRS=\
  MAPMATCHING\
  ID\
  SVM

default:	dirs $(OBJECTS)

dirs:	dummy
	@$(MAKE) subdirs SUBDIRS="$(SECOND_LEVEL_DIRS)" TARGET="default"

subdirs:	dummy
	@set -e; for i in $(SUBDIRS); do cd $$i; if $(MAKE) $(TARGET); then cd .. ; else cd .. ; exit 1 ; fi ; done

dummy:

depend:	dummy
	@echo "" > .Dependencies
	@$(MAKE) subdirs SUBDIRS="$(SECOND_LEVEL_DIRS)" TARGET=depend

default_clean:
	@-$(RM) -rf .Dependencies *~  core.*
	@$(MAKE) subdirs SUBDIRS="$(SECOND_LEVEL_DIRS)" TARGET=default_clean

collectlib:
	@$(MAKE) subdirs SUBDIRS="$(SECOND_LEVEL_DIRS)" TARGET=collectlib


include ../config.mak

