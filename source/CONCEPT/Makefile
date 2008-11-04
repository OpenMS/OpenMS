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

DIRECTORY=CONCEPT


# Note: VersionInfo should be the last target such that ADD_CXXFLAGS will not affect the other targets.
CPP_SOURCES=\
	ClassTest.C\
	Constants.C\
	Exception.C\
	FactoryProduct.C\
	SingletonRegistry.C\
	ProgressLogger.C\
	TimeStamp.C\
	FuzzyStringComparator.C\
	VersionInfo.C\

include ../common.mak
-include .Dependencies


# Always remake VersionInfo.o to reflect the current status of your SVN working copy
VersionInfo.o: always

always:
	@$(eval ADD_CXXFLAGS += -DOPENMS_REVISION=\"$$(shell set LC_ALL=C; cd ../../; svnversion -n)\")
	@$(ECHO) "Note:  VersionInfo.o  is always (re-)compiled to reflect the current status of your SVN working copy"

