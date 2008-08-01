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

all:	default

include ../config_defs.mak

DIRECTORY=COMPARISON


-include .Dependencies

# Take care of subdirectories

default: $(OBJECTS) dummy
	@cd MAPS; if $(MAKE) default; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd SPECTRA; if $(MAKE) default; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd CLUSTERING; if $(MAKE) default; then cd .. ; else cd .. ; exit 1 ; fi ;

depend:	dummy
	@echo "" > .Dependencies
	@cd MAPS; if $(MAKE) depend; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd SPECTRA; if $(MAKE) depend; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd CLUSTERING; if $(MAKE) depend; then cd .. ; else cd .. ; exit 1 ; fi ;

default_clean: dummy
	@-$(RM) -rf .Dependencies *~  core.*
	@cd MAPS; if $(MAKE) default_clean; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd SPECTRA; if $(MAKE) default_clean; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd CLUSTERING; if $(MAKE) default_clean; then cd .. ; else cd .. ; exit 1 ; fi ;

collectlib: dummy
	@cd MAPS; if $(MAKE) collectlib; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd SPECTRA; if $(MAKE) collectlib; then cd .. ; else cd .. ; exit 1 ; fi ;
	@cd CLUSTERING; if $(MAKE) collectlib; then cd .. ; else cd .. ; exit 1 ; fi ;

dummy:

include ../config.mak

