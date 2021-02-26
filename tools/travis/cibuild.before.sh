#!/bin/bash
export CI_PROVIDER="Travis-CI"
if [ "$TRAVIS_OS_NAME" = "linux" ]; then
  bash tools/travis/lnx-cibuild.before.sh
else 
  bash tools/travis/osx-cibuild.before.sh
fi

