# -*- mode: C++; tab-width: 2; -*-
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
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

#Make sure a SVN path is given as first argument
#Make sure the path to the documentation is given as second argument
if test ! $1 || test ! $2; then 
	echo "Use: make_dist.sh <branch> <documentation>";
 	echo "";
 	echo "Pass the SVN path inside the OpenMS repository as first argument.";
 	echo "For the 1.2 release branch that would be 'branches/Release1.2'.";
 	echo "";
 	echo "Pass the path to an OpemMS/doc directory containing the current";
 	echo "documentation as the second argument.";
 	echo "";
 	exit;
fi


# create a dummy directory and extract the branch of OpenMS
echo "####################################################"
echo "extracting OpenMS"
cd /tmp
rm -rf OpenMS-dist
mkdir OpenMS-dist
cd /tmp/OpenMS-dist
svn export https://open-ms.svn.sourceforge.net/svnroot/open-ms/$1 || ( echo "could not extract OpenMS" >&0 && exit )
mv `basename $1` OpenMS
echo ""

# copy documentation
echo "####################################################"
echo "copying documentation to release"
cp $2/OpenMS_tutorial.pdf OpenMS/doc/
cp $2/TOPP_tutorial.pdf OpenMS/doc/
cp -R $2/html OpenMS/doc/

# extract the current SVN version of contrib
echo "####################################################"
echo "extracting contrib"
svn export https://open-ms.svn.sourceforge.net/svnroot/open-ms/contrib || ( echo "cannot extract contrib" >&0 && exit )

# copy contrib
echo "####################################################"
echo "copying contrib to release"
cp -R contrib OpenMS

# extracting version number
echo "####################################################"
echo "extracting version number"
VERSION="`grep AC_INIT /tmp/OpenMS-dist/OpenMS/source/configure.ac | awk -F, '{print 'bla'$2'bla'}' | sed -e 's/\s//g'`"
echo "Version: ${VERSION}"
echo ""

# create archive
echo "####################################################"
echo "creating archive"
DIR="OpenMS-${VERSION}"
FILE="${DIR}.tar.gz"
cd /tmp/OpenMS-dist/
mv OpenMS ${DIR}
echo "creating archive /tmp/OpenMS-dist/${FILE}"
tar zcf ${FILE} ${DIR}
