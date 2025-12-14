#!/usr/bin/env bash

set -eu
set -o pipefail

# Additional packages that we use when building deb packages always runs after deps-ubuntu.sh

sudo apt-get -qq install -y \
   libsqlitecpp-dev\
   nlohmann-json3-dev\
   libsimde-dev