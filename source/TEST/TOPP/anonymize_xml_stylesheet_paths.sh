# -*- mode: shell-script; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
# $Maintainer: Clemens Groepl $
# $Authors: $
# --------------------------------------------------------------------------

# PURPOSE: This will change the path to xml stylesheets.  Thus you can see if
# the whitelist for FuzzyDiff is set such that your TOPP test would pass on
# another platform and for other users.  Apply this to all template files
# which contain paths to xml stylesheets.

# USAGE: anonymize_xml_stylesheet_paths file1 file2 file3


sed -i.sed_backup.`date +%Y%m%d%H%M%S` 's:xml-stylesheet\ type=\"text/xsl\"\ href=\"file\://.*/\(.*\)\":xml-stylesheet\ type=\"text/xsl\"\ href=\"file\:///some/path/to/\1:' $*
