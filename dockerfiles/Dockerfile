# syntax=docker/dockerfile:1.3-labs
# default OpenMS version is "latest" for the current tip of master
# when this image is built from releases it will be "3.1.0" or similar
ARG OPENMS_VERSION_TAG=latest
ARG SOURCE_DIR="/tmp/OpenMS"
ARG BUILD_DIR="${SOURCE_DIR}/bld"
ARG INSTALL_DIR="/opt/OpenMS"
# the contrib image is the collection of dependencies needed to build OpenMS
ARG CONTRIB_IMAGE_TAG=latest
ARG CONTRIB_BASE_IMAGE=ghcr.io/openms/contrib:${CONTRIB_IMAGE_TAG}
# this directory is where any dependencies that are built from source go
ARG CONTRIB_BUILD_DIR="/contrib-build"
ARG THIRDPARTY_DIR="/opt/OpenMS/thirdparty"
ARG NUM_BUILD_CORES=8
ARG MAKEFLAGS="-j${NUM_BUILD_CORES}"
ARG CMAKE_VERSION="3.28.1"


### OpenMS base ###
FROM ${CONTRIB_BASE_IMAGE} AS base
ARG OPENMS_VERSION_TAG
ARG CONTRIB_BUILD_DIR

# if the contrib image builds anything from source it will place the binaries
# in the CONTRIB_BUILD_DIR. If it doesn't (ie. because they are installed from
# the package manager) then we need to create the directory anyway since other
# stages attempt to copy from it
RUN mkdir -p ${CONTRIB_BUILD_DIR}/bin

# the main reason we're pulling this in is to make sure that the metadata
# below is consistently formatted for releases irrespective of the base image
LABEL version="${OPENMS_VERSION_TAG}"
LABEL software.version="${OPENMS_VERSION_TAG}-Ubuntu22.04"
LABEL org.opencontainers.image.source https://github.com/OpenMS/OpenMS


### OpenMS build ###
FROM base AS build
ARG SOURCE_DIR
ARG BUILD_DIR
ARG INSTALL_DIR
ARG CONTRIB_BUILD_DIR
ARG THIRDPARTY_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

COPY . ${SOURCE_DIR}

WORKDIR ${SOURCE_DIR}
RUN mkdir -p ${THIRDPARTY_DIR}
RUN <<-EOF
    set -eux
    git submodule update --depth=1 --init THIRDPARTY

    # copying only the binaries that are relevant to Linux
    cp -r THIRDPARTY/All/* ${THIRDPARTY_DIR}
    cp -r THIRDPARTY/Linux/64bit/* ${THIRDPARTY_DIR}

    rm -rf THIRDPARTY
EOF

# adding thirdparty binaries to the PATH so they will be picked up by the build system
ENV PATH="${THIRDPARTY_DIR}/LuciPHOr2:${THIRDPARTY_DIR}/MSGFPlus:${THIRDPARTY_DIR}/ThermoRawFileParser:${THIRDPARTY_DIR}/Comet:${THIRDPARTY_DIR}/Fido:${THIRDPARTY_DIR}/MaRaCluster:${THIRDPARTY_DIR}/Percolator:${THIRDPARTY_DIR}/SpectraST:${THIRDPARTY_DIR}/XTandem:${THIRDPARTY_DIR}/Sage:${PATH}"
RUN cmake \
    -DCMAKE_BUILD_TYPE='Release' \
    -DCMAKE_PREFIX_PATH="${CONTRIB_BUILD_DIR}/;/usr/;/usr/local" \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DBOOST_USE_STATIC=OFF \
    -S ${SOURCE_DIR} \
    -B ${BUILD_DIR}

WORKDIR ${BUILD_DIR}
RUN make all
RUN make install/strip


### OpenMS library ###
FROM ubuntu:22.04 AS library
ARG INSTALL_DIR
ARG CONTRIB_BUILD_DIR

RUN apt-get update \
  && apt-get install -y --no-install-recommends --no-install-suggests \
    libqt5opengl5 \
    libsvm3 \
    libzip4 \
    zlib1g \
    libbz2-1.0 \
    libgomp1 \
    libqt5svg5 \
    libxerces-c3.2 \
    coinor-libcoinmp1v5 \
    libboost-date-time1.74-dev  \
    libboost-iostreams1.74-dev \
    libboost-regex1.74-dev \
    libboost-math1.74-dev \
    libboost-random1.74-dev \
  && rm -rf /var/lib/apt/lists/*

COPY --from=build ${CONTRIB_BUILD_DIR}/bin ${CONTRIB_BUILD_DIR}/bin
COPY --from=build ${INSTALL_DIR}/lib ${INSTALL_DIR}/lib
COPY --from=build ${INSTALL_DIR}/include ${INSTALL_DIR}/include
COPY --from=build ${INSTALL_DIR}/share ${INSTALL_DIR}/share

LABEL software="OpenMS (library)"


### OpenMS tools runtime ###
FROM library AS tools
ARG INSTALL_DIR

COPY --from=build ${INSTALL_DIR}/bin ${INSTALL_DIR}/bin
ENV PATH="${INSTALL_DIR}/bin:${PATH}"

LABEL software="OpenMS (tools)"


### OpenMS tools runtime ###
FROM tools AS tools-thirdparty
ARG THIRDPARTY_DIR

RUN apt-get update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends --no-install-suggests \
    ## thirdparty dependencies
    # LuciPHOr2 and MSGFPlus:
    openjdk-17-jre-headless \
    # ThermoRawFileParser:
    mono-runtime \
    # NOTE: we need mono-devel because ThermoRawFileParser needs
    #       netstandard.dll, which is not included in mono-runtime.
    #       See: https://github.com/mono/mono/issues/14902,
    #            https://github.com/mono/mono/issues/17148
    #       It's possible that we could re-build it to avoid mono altogther.
    mono-devel \
  && rm -rf /var/lib/apt/lists/*

COPY --from=build ${THIRDPARTY_DIR} ${THIRDPARTY_DIR}

ENV PATH="${THIRDPARTY_DIR}/LuciPHOr2:${THIRDPARTY_DIR}/MSGFPlus:${THIRDPARTY_DIR}/ThermoRawFileParser:${THIRDPARTY_DIR}/Comet:${THIRDPARTY_DIR}/Fido:${THIRDPARTY_DIR}/MaRaCluster:${THIRDPARTY_DIR}/Percolator:${THIRDPARTY_DIR}/SpectraST:${THIRDPARTY_DIR}/XTandem:${THIRDPARTY_DIR}/Sage:${PATH}"

LABEL software="OpenMS (tools + thirdparty)"


### OpenMS runtime test ###
# this stage is based on the runtime stage, so it's a test to ensure
# that the runtime dependencies satisfy the requirements of the tools,
# binaries and the OpenMS library
FROM tools-thirdparty AS test
ARG SOURCE_DIR
ARG BUILD_DIR
ARG CMAKE_VERSION
ARG NUM_BUILD_CORES

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends --no-install-suggests \
    # we need Xvfb to run a small subset of tests (eg. TOPP_INIUpdater)
    xvfb \
    # needed for TSGDialog_test and TOPPView_test
    libqt5test5 \
  && rm -rf /var/lib/apt/lists/*

# installing cmake
WORKDIR /tmp
ADD https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.sh cmake.sh
RUN <<-EOF
    set -eux
    mkdir -p /opt/cmake
    sh cmake.sh --skip-license --prefix=/opt/cmake
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
    ln -s /opt/cmake/bin/ctest /usr/local/bin/ctest
    rm -rf /tmp/*
EOF

COPY --from=build ${SOURCE_DIR} ${SOURCE_DIR}
COPY --from=build ${BUILD_DIR} ${BUILD_DIR}

WORKDIR ${BUILD_DIR}
RUN xvfb-run -a ctest --output-on-failure -j${NUM_BUILD_CORES}
