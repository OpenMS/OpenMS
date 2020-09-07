# Create manylinux packages from current release using the docker image
#
# based on https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
#
# Execute as:
# 
#   sudo docker run --net=host -v `pwd`:/data hroest/manylinux2014_qt59_contrib:v1.1 /bin/bash /data/create-manylinux.sh
#

set -e

# Note that we cannot use git any more as it refuses to communicate with
# github, but we have a patched version of wget with a newer OpenSSL version
# capable of downloading the release tar file.
wget https://github.com/OpenMS/OpenMS/releases/download/Release2.6.0/OpenMS-2.6.0-src.tar.gz -O OpenMS-2.6.0-src.tar.gz
tar xzvf OpenMS-2.6.0-src.tar.gz
mv OpenMS-2.6.0/ OpenMS

# Apply a patches / fixes
# 1. Updated readme (due to pypi.org not displaying markdown)
# 2. Updated version number
# 3. Patch init (since the linked libraries are all in .lib)
cd OpenMS
wget https://raw.githubusercontent.com/hroest/OpenMS/0f3a2e8833d18fc54b7de86e6d67ab7349e828a0/src/pyOpenMS/README.rst -O src/pyOpenMS/README.rst
sed -i 's/@CF_OPENMS_PACKAGE_VERSION@/2.6.0/' src/pyOpenMS/env.py.in
patch -p0 src/pyOpenMS/pyopenms/__init__.py < /data/manylinux.patch 
patch -p0 src/pyOpenMS/setup.py < /data/manylinux.patch2 ## somehow required to link against regex boost lib
cd /

# fixes some issues with empty pxd files
rm -rf /OpenMS/src/pyOpenMS/pxds/SwathMapMassCorrection.pxd


# install Python deps
for PYBIN in /opt/python/cp3*; do
  "$PYBIN/bin/pip" install -U Cython
  "$PYBIN/bin/pip" install -U setuptools
  "$PYBIN/bin/pip" install -U wheel==0.31.1
  "$PYBIN/bin/pip" install -U numpy
  "$PYBIN/bin/pip" install -U nose
  "$PYBIN/bin/pip" install -U autowrap==0.18.1
done

mkdir -p /data/wheelhouse/
mkdir -p /data/wheelhouse/before_fix/

LD_OLD_LIBRARY_PATH=$LD_LIBRARY_PATH

# compile and configure OpenMS
for PYBIN in /opt/python/cp3*; do

  PYVER=`basename $PYBIN`
  mkdir /openms-build-$PYVER
  cd /openms-build-$PYVER
  # configure
  cmake -DCMAKE_PREFIX_PATH="/contrib-build/" -DPYOPENMS=On -DPYTHON_EXECUTABLE:FILEPATH=$PYBIN/bin/python /OpenMS
  make -j4 pyopenms

  # ensure auditwheel can find the libraries
  export LD_LIBRARY_PATH=$LD_OLD_LIBRARY_PATH:`pwd`/lib
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

  /usr/bin/cp wheelhouse/* /data/wheelhouse/before_fix/

  # Somehow auditwheelis broken and we get the wrong libraries here for
  # libcrypto and libssl, lets go and fix it ...
  yum install -y zip
  cd wheelhouse
  mkdir fix
  cp pyopenms-*.whl fix/fix.zip
  cd fix
  unzip fix.zip 
  /usr/bin/cp /lib64/libssl.so.10  pyopenms/.libs/libssl-*
  /usr/bin/cp /lib64/libcrypto.so.10  pyopenms/.libs/libcrypto-*
  zip -r ../pyopenms-*.whl pyopenms/
  cd ../..
  rm -rf wheelhouse/fix/

  mv wheelhouse/* /data/wheelhouse/
done


