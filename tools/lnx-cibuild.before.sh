#!/bin/bash

# fetch contrib and build seqan+gsl
git clone https://github.com/OpenMS/contrib.git
cd contrib
mkdir _build
cd _build

# build seqan from contrib
cmake .. -DBUILD_TYPE=SEQAN
# we build the gsl as the one installed with this ubuntu version
# conflicts with OpenMS
cmake .. -DBUILD_TYPE=GSL -DNUMBER_OF_JOBS=4

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

# build custom cppcheck if we want to perform style tests
if [ $ENABLE_STYLE_TESTING == "On" ]; then
  sudo apt-get install -qq libpcre3-dev
  git clone git://github.com/danmar/cppcheck.git
  pushd cppcheck
  git checkout 1.65
  CXX=clang++ make SRCDIR=build CFGDIR=`pwd`/cfg HAVE_RULES=yes -j4
  popd
else
  # regular builds .. get the search engine executables
  svn checkout http://svn.code.sf.net/p/open-ms/code/THIRDPARTY/SEARCHENGINES/Linux/64bit/ _searchengines
  # remove .svn otherwise we can't check out the other search engines into the same directory (TODO: maybe switch to wget)
  rm _searchengines/.svn -R -f
  svn checkout http://svn.code.sf.net/p/open-ms/code/THIRDPARTY/SEARCHENGINES/All/ _searchengines
fi
