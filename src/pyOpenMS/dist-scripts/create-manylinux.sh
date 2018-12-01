# Create manylinux packages from current release using the docker image
#
# based on https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
#
# Execute as:
# 
#   sudo docker run --net=host -v `pwd`:/data hroest/manylinux_qt59_contrib:v3 /bin/bash /data/create-manylinux.sh
#

set -e

yum -y install zip

git clone -b Release2.4.0 https://github.com/OpenMS/OpenMS.git

# Apply a patches / fixes
# 1. Update version number
# 2. Patch setup.py to ensure we use single-threaded version (3.6/3.7 stalling bug)
# 3. Patch init (since the linked libraries are all in .lib)
cd OpenMS
sed -i 's/@CF_OPENMS_PACKAGE_VERSION@/2.4.0/' src/pyOpenMS/env.py.in
sed -i 's/single_threaded = False/single_threaded = True/' src/pyOpenMS/setup.py
patch -p0 src/pyOpenMS/pyopenms/__init__.py < /data/manylinux.patch 
cd /

# install Python deps
for PYBIN in /opt/python/cp27* /opt/python/cp3[4-9]*; do
  "$PYBIN/bin/pip" install -U Cython
  "$PYBIN/bin/pip" install -U setuptools
  "$PYBIN/bin/pip" install wheel==0.31.1 # bug in newer wheel
  "$PYBIN/bin/pip" install -U numpy
  "$PYBIN/bin/pip" install -U nose
  "$PYBIN/bin/pip" install -U autowrap
done

mkdir -p /data/wheelhouse/

# compile and configure OpenMS
for PYBIN in /opt/python/cp27* /opt/python/cp3[4-9]*; do

  PYVER=`basename $PYBIN`
  mkdir /openms-build-$PYVER
  cd /openms-build-$PYVER
  # configure (no OpenMP as this can lead to issues with Python)
  cmake -DCMAKE_PREFIX_PATH="/contrib-build/;/qt/" -DPYOPENMS=On -DPYTHON_EXECUTABLE:FILEPATH=$PYBIN/bin/python \
    -DMT_ENABLE_OPENMP=OFF -DPY_NO_OUTPUT=Off -DPY_SINGLE_THREADED=On \
    -DQT_QMAKE_EXECUTABLE=/qt/bin/qmake ../OpenMS
  make -j6 pyopenms

  # ensure auditwheel can find the libraries
  export LD_LIBRARY_PATH=`pwd`/lib

  # strip the libraries
  strip --strip-all lib/libOpenMS.so 
  strip --strip-all lib/libOpenSwathAlgo.so 
  strip --strip-all lib/libSuperHirn.so 
  cd pyOpenMS
  # remove the libraries as auditwheel will take care of linked libs
  rm -rf pyopenms/lib*
  rm -rf build/lib*/pyopenms/lib*
  "$PYBIN/bin/pip" wheel . -w wheelhouse_tmp

  # Bundle external shared libraries into the wheels
  for whl in wheelhouse_tmp/pyopenms*.whl; do
    auditwheel repair "$whl" -w wheelhouse/
  done

  # We need to fix the Qt libraries as auditwheel seems to strip the qt
  # libraries which produces problems
  unalias cp
  for WHEEL in wheelhouse/*.whl; do
    mkdir -p wheelhouse/fixed
    mkdir /tmp/pytmp
    cp $WHEEL /tmp/pytmp/
    WHEEL=`basename $WHEEL`
    cp /qt/lib/libQt5Core.so /tmp/pytmp/
    cp /qt/lib/libQt5Network.so /tmp/pytmp/
    cd /tmp/pytmp
    unzip $WHEEL
    # Replace the Qt libraries with the original ones
    cp libQt5Core.so pyopenms/.libs/libQt5Core*
    cp libQt5Network.so pyopenms/.libs/libQt5Network*
    rm -rf pyopenms/lib*.so
    zip -r $WHEEL pyopenms/ pyopenms*dist*/
    cd -
    cp /tmp/pytmp/$WHEEL wheelhouse/fixed/
    rm -rf /tmp/pytmp
  done

  mv wheelhouse/fixed/* /data/wheelhouse/
done

