#!/bin/sh

# simple script to configure a build environment for coverity
mkdir _coverity_build
cd _coverity_build

cmake -D CMAKE_FIND_ROOT_PATH="${PWD}/contrib/_build/;/usr;" -D BOOST_USE_STATIC=Off -D CMAKE_BUILD_TYPE=Release ..
