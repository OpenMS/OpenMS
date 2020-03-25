#!/bin/bash
# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry               
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
# $Maintainer: Julianus Pfeuffer$
# $Authors: Marc Sturm, Julianus Pfeuffer $
# --------------------------------------------------------------------------

# Takes source and contrib (freshly checked out from given branch or by specifying git folders),
# as well as a prebuilt documentation (should include tutorials), copies it
# to one folder and 

while [[ $# > 1 ]]
do
key="$1"

case $key in
    -s|--source)
    SOURCEPATH="$2"
    shift # past argument
    ;;
    -c|--contrib)
    CONTRIBPATH="$2"
    shift # past argument
    ;;
    -d|--docs)
    DOCPATH="$2"
    shift # past argument
    ;;
    -b|--branch)
    BRANCH="$2"
    shift # past argument
    ;;
    -r|--releaseversion)
    RELVERSION="$2"
    shift # past argument
    ;;
    *)
    # unknown option
    ;;
esac

shift # past argument or value
done

if [[ -z $BRANCH ]]; then
    # branch defaults to master
    BRANCH=master
fi

if [[ -z $CONTRIBPATH ]]; then
    # Following line: to work on LNX and OSX
    mycontribtmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mycontribtmpdir'`

    echo "Cloning contrib master to $mycontribtmpdir/contrib"
    # needs git
    git clone https://github.com/OpenMS/contrib/ $mycontribtmpdir/contrib || {
     echo "Can not clone contrib to $mycontribtmpdir/contrib. Specify an already cloned contrib path via -c or make sure the folder is empty. Also check internet connection";
     exit ${LINENO} 
    } >&2
    
    CONTRIBPATH="$mycontribtmpdir/contrib"
fi

# If source not given, clone here
if [[ -z $SOURCEPATH ]]; then
    # Following line: to work on LNX and OSX
    mysourcetmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mysourcetmpdir'`

    echo "Cloning OpenMS sources of branch $BRANCH to $mysourcetmpdir/OpenMS-$RELVERSION"

    # needs git
    git clone https://github.com/OpenMS/OpenMS/ $mysourcetmpdir/OpenMS-$RELVERSION || {
     echo "Can not clone sources to $mysourcetmpdir/OpenMS-$RELVERSION. Specify an already cloned OpenMS source path via -s or make sure the folder is empty. Also check internet connection";
     exit ${LINENO} 
    } >&2
    
    SOURCEPATH="$mysourcetmpdir/OpenMS-$RELVERSION"
fi

if [[ -z $DOCPATH ]]; then
     echo "Path to documentation not specified. Please build the documentation with tutorials and specify path with -d";
     exit ${LINENO}
fi

## Release version can be empty. Only affects name of the tarball

echo Used the following variables:
echo SOURCE  = "${SOURCEPATH}"
echo CONTRIB     = "${CONTRIBPATH}"
echo DOCU    = "${DOCPATH}"
echo BRANCH    = "${BRANCH}"
echo RELVER  = "${RELVERSION}"

# copy documentation
echo "####################################################"
echo "copying documentation to release"
cp $DOCPATH/OpenMS_tutorial.pdf $SOURCEPATH/doc/
cp $DOCPATH/TOPP_tutorial.pdf $SOURCEPATH/doc/
cp -r $DOCPATH/html $SOURCEPATH/doc/

# move contrib
echo "####################################################"
echo "moving contrib into sources"
cp -r $CONTRIBPATH $SOURCEPATH

# create archive
echo "####################################################"
echo "creating archive"
TARNAME="OpenMS-$RELVERSION"
FILE="${TARNAME}-sources.tar.gz"
echo "creating archive $PWD/${FILE}"
tar zcf ${FILE} -C ${SOURCEPATH} ..
