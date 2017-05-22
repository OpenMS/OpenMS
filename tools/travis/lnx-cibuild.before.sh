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

# necessary for succesfull build?
build_contrib COINOR

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
  svn checkout https://github.com/OpenMS/THIRDPARTY/trunk/Linux/64bit/ _thirdparty
  # remove .svn otherwise we can't check out the other search engines into the same directory
  rm _thirdparty/.svn -R -f || true
  svn checkout https://github.com/OpenMS/THIRDPARTY/trunk/All/ _thirdparty
fi


