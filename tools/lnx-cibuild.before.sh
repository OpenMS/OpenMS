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
cd contrib

# we build seqan as the versions shipped in Ubuntu are not recent enough
build_contrib SEQAN

# we build the gsl as the one installed with this ubuntu version
# conflicts with OpenMS
build_contrib GSL

# add alternative repo for newer boost version
sudo add-apt-repository --yes ppa:boost-latest/ppa
sudo apt-get update

# install required packages
sudo apt-get install -qq  boost1.54\
                          libxerces-c3.1\
                          libxerces-c-dev \
                          libicu-dev \
                          qt4-dev-tools \
                          libqt4-dev \
                          libqt4-core \
                          libqt4-gui \
                          libsvm-dev \
                          libsvm3 \
                          glpk

# install doxygen but not latex (saves some time)
sudo apt-get install -qq \
                     --no-install-recommends doxygen \
                     graphviz
