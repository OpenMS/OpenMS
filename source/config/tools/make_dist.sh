# create a dummy directory and extract the current SVN version of OpenMS
echo "extracting OpenMS"
cd /tmp
rm -rf OpenMS-dist
mkdir OpenMS-dist
cd /tmp/OpenMS-dist
svn co https://svn.sourceforge.net/svnroot/open-ms/OpenMS 1>svn_extract.log 2>svn_extract_err.log || ( echo "cannot extract OpenMS" >&0 && exit )

# remove SVN information
echo "removing SVN info"
REMOVE=`find . -name .svn -type d`
rm -rf ${REMOVE} 2>/dev/null

# create archive
DIR="OpenMS-`grep AC_INIT /tmp/OpenMS-dist/OpenMS/source/configure.ac | awk -F, '{print 'bla'$2'bla'}' | sed -e 's/\s//g'`"
FILE="${DIR}.tar.gz"
cd /tmp/OpenMS-dist/
mv OpenMS ${DIR}

echo "creating archive /tmp/OpenMS-dist/${FILE}"
tar zcf ${FILE} ${DIR}
