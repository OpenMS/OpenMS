# Create manylinux packages from current release using the docker image
#
# based on https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
#
# Execute as:
# 
#   sudo docker run --net=host -v `pwd`:/data hroest/manylinux_qt_contrib:v2 /bin/bash /data/create-manylinux.sh
#

set -e

# Note that we cannot use git any more as it refuses to communicate with
# github, but we have a patched version of wget with a newer OpenSSL version
# capable of downloading the release tar file.
wget https://github.com/OpenMS/OpenMS/releases/download/Release2.3.0/OpenMS-2.3.0-src.tar.gz -O OpenMS-2.3.0-src.tar.gz
tar xzvf OpenMS-2.3.0-src.tar.gz
mv OpenMS-2.3.0/ OpenMS

# Apply a patches / fixes
# 1. Updated readme (due to pypi.org not displaying markdown)
# 2. Updated version number
# 3. Patch init (since the linked libraries are all in .lib)
cd OpenMS
wget https://raw.githubusercontent.com/hroest/OpenMS/0f3a2e8833d18fc54b7de86e6d67ab7349e828a0/src/pyOpenMS/README.rst -O src/pyOpenMS/README.rst
sed -i 's/@CF_OPENMS_PACKAGE_VERSION@/2.3.0.4/' src/pyOpenMS/env.py.in
patch -p0 src/pyOpenMS/pyopenms/__init__.py < /data/manylinux.patch 
cd /


# install Python deps
for PYBIN in /opt/python/cp27* /opt/python/cp3[4-9]*; do
  "$PYBIN/bin/pip" install -U Cython
  "$PYBIN/bin/pip" install -U setuptools
  "$PYBIN/bin/pip" install -U wheel
  "$PYBIN/bin/pip" install -U numpy
  "$PYBIN/bin/pip" install -U nose
  "$PYBIN/bin/pip" install -U autowrap
done

mkdir /data/wheelhouse/

# compile and configure OpenMS
for PYBIN in /opt/python/cp27* /opt/python/cp3[4-9]*; do

  PYVER=`basename $PYBIN`
  mkdir /openms-build-$PYVER
  cd /openms-build-$PYVER
  # configure
  cmake -DCMAKE_PREFIX_PATH="/contrib-build/" -DPYOPENMS=On -DPYTHON_EXECUTABLE:FILEPATH=$PYBIN/bin/python -DQT_QMAKE_EXECUTABLE=/qt/bin/qmake ../OpenMS
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
  mv wheelhouse/* /data/wheelhouse/
done


