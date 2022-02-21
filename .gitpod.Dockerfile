FROM ghcr.io/openms/contrib:latest
                    
#USER gitpod
# Avoid user input
ARG DEBIAN_FRONTEND=noninteractive
ARG CLANGDVER=12.0.1
ARG CMAKE_VERSION=3.22.2

RUN wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh \
      && chmod u+x /tmp/cmake-install.sh \
      && mkdir /usr/bin/cmake \
      && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
      && rm /tmp/cmake-install.sh

ENV PATH="/usr/bin/cmake/bin:${PATH}"

# Install tools for VSCode, Intellisense, JRE for Thirdparties, etc. using apt-get
RUN apt-get -q update && apt-get install -yq gdb unzip wget php openjdk-11-jre python3-pip && rm -rf /var/lib/apt/lists/*
RUN wget https://github.com/clangd/clangd/releases/download/$CLANGDVER/clangd-linux-$CLANGDVER.zip && unzip clangd-linux-$CLANGDVER.zip -d /opt/ && mv /opt/clangd_$CLANGDVER /opt/clangd/ && rm clangd-linux-$CLANGDVER.zip 
#
# More information: https://www.gitpod.io/docs/42_config_docker/
