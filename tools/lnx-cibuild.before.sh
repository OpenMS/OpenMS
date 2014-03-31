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



# leave contrib
popd

# add alternative repo for newer boost version
sudo add-apt-repository --yes ppa:boost-latest/ppa
sudo apt-get update

# install required packages
sudo apt-get install -qq  libboost-date-time1.54-dev \
                          libboost-iostreams1.54-dev \
                          libboost-regex1.54-dev \
                          libboost-math1.54-dev \
                          libboost-random1.54-dev \
                          libxerces-c3.1 \
                          libxerces-c-dev \
                          libicu-dev \
                          qt4-dev-tools \
                          libqt4-dev \
                          libqt4-core \
                          libqt4-gui \
                          libsvm-dev \
                          libsvm3 \
                          glpk \
                          subversion                             

# install doxygen but not latex (saves some time)
sudo apt-get install -qq \
                     --no-install-recommends doxygen \
                     graphviz

# get the search engine executables 
svn checkout http://svn.code.sf.net/p/open-ms/code/THIRDPARTY/SEARCHENGINES/Linux/64bit/ _searchengines
