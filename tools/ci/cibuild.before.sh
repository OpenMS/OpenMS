#!/bin/bash
if [ "$TRAVIS_OS_NAME" = "linux" ]; then
  bash tools/travis/lnx-cibuild.before.sh
else 
  bash tools/travis/osx-cibuild.before.sh
fi

