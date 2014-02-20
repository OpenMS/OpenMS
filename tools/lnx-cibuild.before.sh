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

# packages to install

# install required packages
# we only need a subset for coverity, since we want to 
# documentation there
if [ "$COVERITY_SCAN_BRANCH" == "1" ]; then
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
else
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
	                          glpk \
	                          doxygen \
	                          graphviz
fi
