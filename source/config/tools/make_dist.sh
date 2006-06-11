#  $Id: make_dist.sh,v 1.2 2006/05/26 09:07:50 marc_sturm Exp $

# create a dummy directory and extract the current SVN version of OpenMS
echo "extracting OpenMS ..."
cd /tmp
rm -rf OpenMS-dist
mkdir OpenMS-dist
cd OpenMS-dist
svn co http://svn.sourfeforge.net/svnroot/open-ms/OpenMS 1>svn_extract.log 2>svn_extract_err.log || ( echo "cannot extract OpenMS" >&0 && exit )
# remove SVN information
REMOVE=`find . -name .svn -type d`
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
