#!/usr/bin/env bash

# Create manylinux packages from current release using the docker image
#
# based on https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh
#
# Execute as:
#
#   sudo docker run --net=host -v `pwd`:/data hroest/manylinux_qt59_contrib:v1.2 /bin/bash /data/create-manylinux.sh
#

## For a release, change to the following:
## git clone -b Release2.3.0 https://github.com/OpenMS/OpenMS.git
git clone https://github.com/OpenMS/OpenMS.git 
 
# Bugfix 1:
# make sure that we can find the link library
ln -s /contrib-build/lib64/libxerces-c-3.2.a /contrib-build/lib/libxerces-c.a

yum install -y zip

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

  PYVER=$(basename $PYBIN)
  mkdir /openms-build-$PYVER
  cd /openms-build-$PYVER

  # configure and build
  cmake -DCMAKE_PREFIX_PATH="/qt/;/contrib-build/" -DPYOPENMS=On -DPYTHON_EXECUTABLE:FILEPATH=$PYBIN/bin/python \
    -DPY_NUM_THREADS=2 -DPY_NUM_MODULES=4 -DQT_QMAKE_EXECUTABLE=/qt/bin/qmake ../OpenMS
  make -j6 pyopenms

  # create final wheel ready for bundling (and make sure we find the OpenMS libs)
  export LD_LIBRARY_PATH=$(pwd)/lib
  cd pyOpenMS
  "$PYBIN/bin/pip" wheel . -w wheelhouse_tmp

  # Bundle external shared libraries into the wheels
  for whl in wheelhouse_tmp/pyopenms*.whl; do
    auditwheel repair "$whl" -w wheelhouse/
  done

  # Changing the soname of libQt5Core does not end well and leads to segfaults
  # when using the module.
  # Therefore we copy the original libQt5Core back into the wheel
  mkdir wheelhouse_fixed
  for whl in wheelhouse/pyopenms*.whl; do
    /bin/cp $whl wheelhouse_fixed
    whln=`basename $whl`
    bn=`basename $whl .whl`
    cd wheelhouse_fixed
    mv $bn.whl $bn.zip
    unzip $bn.zip
    rm -rf $bn.zip
    /bin/cp /qt/lib/libQt5Core.so pyopenms/.libs/libQt5Core-*
    rm -rf pyopenms/lib*
    # this does not always work, see https://github.com/pypa/manylinux/issues/119
    # strip --strip-unneeded pyopenms/.libs/libOpenMS-*.so
    zip -r $bn.zip pyopenms pyopenms-*.dist-info
    rm -rf pyopenms pyopenms-*.dist-info
    mv $bn.zip $bn.whl
    cd ..
  done

  # retrieve data
  mv wheelhouse_fixed/* /data/wheelhouse/

done

