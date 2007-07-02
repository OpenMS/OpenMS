# -*- Mode: C++; tab-width: 2; -*-
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

#Make sure a SVN path is given
 if test ! $1; then 
 	echo "Please pass the SVN path inside the OpenMS repository as first argument!";
 	echo "For the HEAD revision that would be 'OpenMS'.";
 	exit;
 fi

# create a dummy directory and extract the current SVN version of OpenMS
echo "extracting OpenMS"
cd /tmp
rm -rf OpenMS-dist
mkdir OpenMS-dist
cd /tmp/OpenMS-dist
svn co https://svn.sourceforge.net/svnroot/open-ms/$1 1>svn_extract.log 2>svn_extract_err.log || ( echo "cannot extract OpenMS" >&0 && exit )
mv `basename $1` OpenMS
echo ""

# extracting version number
echo "extracting version number"
VERSION="`grep AC_INIT /tmp/OpenMS-dist/OpenMS/source/configure.ac | awk -F, '{print 'bla'$2'bla'}' | sed -e 's/\s//g'`"
echo "Version: ${VERSION}"
echo ""

# remove SVN information
echo "removing SVN info"
REMOVE=`find . -name .svn -type d`
rm -rf ${REMOVE} 2>/dev/null
echo ""

# create archive
DIR="OpenMS-${VERSION}"
FILE="${DIR}.tar.gz"
cd /tmp/OpenMS-dist/
mv OpenMS ${DIR}
echo "creating archive /tmp/OpenMS-dist/${FILE}"
tar zcf ${FILE} ${DIR}
