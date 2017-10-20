#!/bin/bash

function build_contrib {
  cmake . -DBUILD_TYPE=$1

  if [ $? -ne 0 ]; then
    # we give it another try
    echo "1st attempt to build $1 failed .. retry"
    cmake . -DBUILD_TYPE=$1 -DNUMBER_OF_JOBS=4

    if [ $? -ne 0 ]; then
      echo "2nd attempt to build $1 failed .. abort"
      exit $?
    fi
  fi
}

if [ "${PYOPENMS}" = "ON" ]; then
  # Note: ensure that cmake uses the same python!
  pyenv versions
  # select the desired Python version
  pyenv global 2.7.13
  which pip
  which python

  # small patch to accelerate build
  pwd
  ls src/pyOpenMS/create_cpp_extension.py
  sed -i 's/import time/import time\nPY_NUM_THREADS=4/g' src/pyOpenMS/create_cpp_extension.py

  pip install -U setuptools
  pip install -U pip
  pip install -U nose
  pip install -U numpy
  pip install -U wheel
  pip install -U Cython
  pip install -U autowrap==0.14
fi

# fetch contrib and build seqan
git clone git://github.com/OpenMS/contrib/
pushd contrib

# we build seqan as the versions shipped in Ubuntu are not recent enough
build_contrib SEQAN

# we build WildMagic
build_contrib WILDMAGIC

# we build Eigen as the versions shipped in Ubuntu are not recent enough
build_contrib EIGEN

# we build Sqlite as the versions shipped in Ubuntu are not recent enough
build_contrib SQLITE

# leave contrib
popd


# build custom cppcheck if we want to perform style tests
if [ "${ENABLE_STYLE_TESTING}" = "ON" ]; then
  git clone git://github.com/danmar/cppcheck.git
  pushd cppcheck
  git checkout 1.65
  CXX=clang++ make SRCDIR=build CFGDIR=`pwd`/cfg HAVE_RULES=yes -j4
  popd
else
  # regular builds .. get the search engine executables via githubs SVN interface (as git doesn't allow single folder checkouts)
  svn export --force https://github.com/OpenMS/THIRDPARTY/trunk/Linux/64bit/ _thirdparty
  svn export --force https://github.com/OpenMS/THIRDPARTY/trunk/All/ _thirdparty
fi


