#  $Id: make_dist.sh,v 1.2 2006/05/26 09:07:50 marc_sturm Exp $

# check whether the revision number is given as the first argument
# (CVS tag!)
if test $# != 2 ; then
	echo "make_dist <CVS tag> <SF login>"
	exit
else
	TAG=$1
	LOGIN=$2
fi

# create a dummy directory and extract the current CVS version of OpenMS
echo "extracting OpenMS revision $1 from cvs..."
cd /tmp
rm -rf OpenMS-dist
mkdir OpenMS-dist
cd OpenMS-dist
cvs -q -d:ext:${LOGIN}@open-ms.cvs.sourceforge.net:/cvsroot/open-ms co -r ${TAG} OpenMS 1>cvs_extract.log 2>cvs_extract_err.log || ( echo "cannot extract revision ${TAG} of OpenMS" >&0 && exit )
# remove CVS information
REMOVE=`find . -name CVS -type d`
rm -rf ${REMOVE} 2>/dev/null
cd OpenMS/source
FILE="OpenMS-`grep OPENMS_RELEASE_STRING ../include/OpenMS/CONCEPT/VersionInfo.h | head -1 | awk '{print $3}'| tr -d \\"`.tar"
cd ../../

echo "checking versions..."
grep "e OPENMS_RELEASE_STRING" OpenMS/include/OpenMS/CONCEPT/VersionInfo.h
grep "PROJECT_NUMBER" OpenMS/doc/Doxyfile
grep "@version" OpenMS/include/OpenMS/OpenMS.doxygen
grep "LIB_VERSION" OpenMS/source/config/Makefile.in
grep "AC_INIT" OpenMS/source/configure.ac

echo "creating archive ${FILE}..."
DIR=`basename ${FILE} .tar` 
mv OpenMS ${DIR}
tar zcf ${FILE} ${DIR}
