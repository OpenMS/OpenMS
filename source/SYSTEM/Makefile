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

include ../config_defs.mak

DIRECTORY=SYSTEM

# Where the headers of the derived classes reside
HEADER_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)

# Where the UIC-created headers reside
UIC_DIR=$(OPENMS_PATH)/include/OpenMS/$(DIRECTORY)/UIC

CPP_SOURCES=\
	FileWatcher.C\
	StopWatch.C\
	File.C\
	ProcessResource.C\
	ExternalAllocator.C\
	ExternalAllocatorUnique.C

MOC_SOURCES=\
	FileWatcher_moc.C

SOURCES		= $(CPP_SOURCES) $(MOC_SOURCES)
OBJECTS 	= $(SOURCES:.C=.o)

-include .Dependencies
include ../common.mak


FileWatcher_moc.C:  $(HEADER_DIR)/FileWatcher.h
	$(MOC) $(HEADER_DIR)/FileWatcher.h -o FileWatcher_moc.C
