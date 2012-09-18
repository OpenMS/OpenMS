# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry               
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
# 
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution 
#    may be used to endorse or promote products derived from this software 
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS. 
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# --------------------------------------------------------------------------
# $Maintainer: $
# $Authors: Marc Sturm $
# --------------------------------------------------------------------------

#Make sure a SVN path is given as first argument
#Make sure the path to the documentation is given as second argument
if test ! $1 || test ! $2 || test ! $3; then 
	echo "Use: make_dist.sh <branch> <documentation> <version>";
 	echo "";
 	echo "Pass the SVN path inside the OpenMS repository as first argument.";
 	echo "For the 1.2 release branch that would be 'branches/Release1.2'.";
 	echo "";
 	echo "Pass the path to an OpenMS/doc directory containing the current";
 	echo "documentation as the second argument.";
 	echo "";
 	echo "Pass the release version as third argument.";
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
mv $2/OpenMS_tutorial.pdf OpenMS/doc/
mv $2/TOPP_tutorial.pdf OpenMS/doc/
mv $2/html OpenMS/doc/

# extract the current SVN version of contrib
echo "####################################################"
echo "extracting contrib"
cd OpenMS
svn export https://open-ms.svn.sourceforge.net/svnroot/open-ms/contrib || ( echo "cannot extract contrib" >&0 && exit )

# create archive
echo "####################################################"
echo "creating archive"
DIR="OpenMS-$3"
FILE="${DIR}.tar.gz"
cd /tmp/OpenMS-dist/
mv OpenMS ${DIR}
echo "creating archive /tmp/OpenMS-dist/${FILE}"
tar zcf ${FILE} ${DIR}
