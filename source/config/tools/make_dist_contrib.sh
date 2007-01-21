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

# create a dummy directory and extract the current SVN version of contrib
echo "extracting contrib"
cd /tmp
rm -rf contrib-dist
mkdir contrib-dist
cd /tmp/contrib-dist
svn co https://svn.sourceforge.net/svnroot/open-ms/contrib 1>svn_extract.log 2>svn_extract_err.log || ( echo "cannot extract contrib" >&0 && exit )

# remove SVN information
echo "removing SVN info"
REMOVE=`find . -name .svn -type d`
rm -rf ${REMOVE} 2>/dev/null

# create archive
DIR="OpenMS_contrib-`grep AC_INIT /tmp/contrib-dist/contrib/configure.ac | awk -F, '{print 'bla'$2'bla'}' | sed -e 's/\s//g'`"
FILE="${DIR}.tar.gz"
cd /tmp/contrib-dist/
mv contrib ${DIR}

echo "creating archive /tmp/contrib-dist/${FILE}"
tar zcf ${FILE} ${DIR}
