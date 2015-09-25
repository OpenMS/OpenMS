#!/bin/bash
if [ "$TRAVIS_OS_NAME" = "linux" ]; then
  bash tools/travis/lnx-cibuild.sh
elif [ "$TRAVIS_OS_NAME" = "osx" ]; then
  bash tools/travis/osx-cibuild.sh
fi

