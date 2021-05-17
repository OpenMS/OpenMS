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


# tap the science tap
brew tap homebrew/science
brew tap homebrew/versions
brew update
brew install libsvm xerces-c glpk coinmp
brew install qt4

# fetch contrib
git clone git://github.com/OpenMS/contrib/
pushd contrib

# we build WildMagic
build_contrib WILDMAGIC

# we build Eigen as the versions shipped in Ubuntu are not recent enough
build_contrib EIGEN

# leave contrib
popd

# regular builds .. get the search engine executables via githubs SVN interface (as git doesn't allow single folder checkouts)
svn export --force http://github.com/OpenMS/THIRDPARTY/trunk/MacOS/64bit/ _thirdparty || true
svn export --force http://github.com/OpenMS/THIRDPARTY/trunk/All/ _thirdparty || true
