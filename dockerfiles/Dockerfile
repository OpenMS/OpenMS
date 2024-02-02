# syntax=docker/dockerfile:1.3-labs
# default OpenMS version is "latest" for the current tip of master
# when this image is built from releases it will be "3.1.0" or similar
ARG OPENMS_VERSION=latest

# this is the directory where the OpenMS git repository gets copied to
ARG OPENMS_SOURCE_DIR="/tmp/OpenMS"

# this is the directory where the OpenMS library and tools are built
ARG OPENMS_LIBRARY_BUILD_DIR="/opt/OpenMS"
ARG OPENMS_LIBRARY_DIR="/opt/OpenMS/lib"

# the contrib image is the collection of dependencies needed to build OpenMS
ARG CONTRIB_IMAGE_TAG=latest
ARG CONTRIB_BASE_IMAGE=ghcr.io/openms/contrib:${CONTRIB_IMAGE_TAG}

# this directory is where any dependencies that are built from source go
ARG CONTRIB_BUILD_DIR="/contrib-build"
ARG THIRDPARTY_BRANCH="master"
ARG THIRDPARTY_DIR="/opt/OpenMS/thirdparty"
ARG NUM_BUILD_CORES=8
ARG MAKEFLAGS="-j${NUM_BUILD_CORES}"
ARG CMAKE_VERSION="3.28.1"


### OpenMS contrib ###
# this is an alternative to the contrib image built in the
# OpenMS/contrib repository. Useful mostly for testing.
FROM ubuntu:22.04 AS contrib
ARG OPENMS_VERSION
ARG CONTRIB_BUILD_DIR
ARG CMAKE_VERSION

# install build dependencies
RUN apt-get -y update \
  && apt-get install -y --no-install-recommends --no-install-suggests \
    # build system dependencies
    g++ \
    make \
    git \
    ca-certificates \
    # OpenMS build dependencies
    libsvm-dev \
    libglpk-dev \
    libzip-dev \
    zlib1g-dev \
    libxerces-c-dev \
    libbz2-dev \
    libomp-dev \
    libhdf5-dev\
    libboost-date-time1.74-dev \
    libboost-iostreams1.74-dev \
    libboost-regex1.74-dev \
    libboost-math1.74-dev \
    libboost-random1.74-dev \
    qtbase5-dev \
    libqt5svg5-dev \
    libqt5opengl5-dev \
    libeigen3-dev \
    coinor-libcoinmp-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates

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

# note that in this image, we're not building any thirdparty dependencies
# as they are all installed from apt. However, we still need to create the
# directory as images that are based on this assume it is there.
RUN mkdir -p ${CONTRIB_BUILD_DIR}/bin
WORKDIR /

LABEL base.image="ubuntu:22.04"
LABEL version="${OPENMS_VERSION}"
LABEL software="OpenMS (dependencies)"
LABEL software.version="${OPENMS_VERSION}-Ubuntu22.04"
LABEL description="Base image to build OpenMS: C++ libraries and tools for MS/MS data analysis"
LABEL website="http://www.openms.org/"
LABEL documentation="http://www.openms.org/"
LABEL license="http://www.openms.org/"
LABEL tags="Proteomics"


### OpenMS base ###
FROM ${CONTRIB_BASE_IMAGE} as base
ARG OPENMS_VERSION

# the main reason we're pulling this in is to make sure that the metadata
# below is consistently formatted for releases irrespective of the base image
LABEL version="${OPENMS_VERSION}"
LABEL software.version="${OPENMS_VERSION}-Ubuntu22.04"
LABEL org.opencontainers.image.source https://github.com/OpenMS/OpenMS


### OpenMS library build ###
FROM base as library-build
ARG OPENMS_SOURCE_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG CONTRIB_BUILD_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

COPY . ${OPENMS_SOURCE_DIR}

RUN cmake \
    -DCMAKE_BUILD_TYPE='Release' \
    -DCMAKE_PREFIX_PATH="${CONTRIB_BUILD_DIR}/;/usr/;/usr/local" \
    -DBOOST_USE_STATIC=OFF \
    -S ${OPENMS_SOURCE_DIR} \
    -B ${OPENMS_LIBRARY_BUILD_DIR}
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN make OpenMS


### OpenMS tools build ###
FROM library-build as tools-build
ARG OPENMS_SOURCE_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG THIRDPARTY_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

WORKDIR ${OPENMS_SOURCE_DIR}
RUN mkdir -p ${THIRDPARTY_DIR}
RUN <<-EOF
    set -eux
    git submodule update --depth=1 --init THIRDPARTY

    # copying only the binaries that are relevant to Linux
    cp -r THIRDPARTY/All/* ${THIRDPARTY_DIR}
    cp -r THIRDPARTY/Linux/64bit/* ${THIRDPARTY_DIR}

    rm -rf THIRDPARTY
EOF

WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN make TOPP


### OpenMS library ###
FROM ubuntu:22.04 AS library
ARG OPENMS_SOURCE_DIR
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
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

COPY --from=library-build ${CONTRIB_BUILD_DIR}/bin ${CONTRIB_BUILD_DIR}/bin
COPY --from=library-build ${OPENMS_LIBRARY_DIR} ${OPENMS_LIBRARY_DIR}
COPY --from=library-build ${OPENMS_SOURCE_DIR}/share ${OPENMS_LIBRARY_BUILD_DIR}/share

LABEL software="OpenMS (library)"


### OpenMS tools runtime ###
FROM library AS tools
ARG OPENMS_LIBRARY_BUILD_DIR

COPY --from=tools-build ${OPENMS_LIBRARY_BUILD_DIR}/bin ${OPENMS_LIBRARY_BUILD_DIR}/bin
ENV PATH="${OPENMS_LIBRARY_BUILD_DIR}/bin:${PATH}"
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}

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

COPY --from=tools-build ${THIRDPARTY_DIR} ${THIRDPARTY_DIR}

ENV PATH="${THIRDPARTY_DIR}/LuciPHOr2:${THIRDPARTY_DIR}/MSGFPlus:${THIRDPARTY_DIR}/ThermoRawFileParser:${THIRDPARTY_DIR}/Comet:${THIRDPARTY_DIR}/Fido:${THIRDPARTY_DIR}/MaRaCluster:${THIRDPARTY_DIR}/Percolator:${THIRDPARTY_DIR}/SpectraST:${THIRDPARTY_DIR}/XTandem:${THIRDPARTY_DIR}/Sage:${PATH}"

LABEL software="OpenMS (tools + thirdparty)"


### OpenMS test build ###
FROM tools-build as test-build
ARG OPENMS_LIBRARY_BUILD_DIR

WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
# TODO: determine if we need to make all or if we can get away with less
RUN make all


### OpenMS runtime test ###
# this stage is based on the runtime stage, so it's a test to ensure
# that the runtime dependencies satisfy the requirements of the tools
# binaries and the OpenMS library
FROM tools AS test
ARG OPENMS_SOURCE_DIR
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG CMAKE_VERSION

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

COPY --from=test-build ${OPENMS_SOURCE_DIR} ${OPENMS_SOURCE_DIR}
COPY --from=test-build ${OPENMS_LIBRARY_BUILD_DIR} ${OPENMS_LIBRARY_BUILD_DIR}
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN <<-EOF
    set -eux
    Xvfb :0 &
    DISPLAY=:0 ctest
EOF


### pyOpenMS ###
FROM library AS pyopenms
ARG OPENMS_LIBRARY_BUILD_DIR
ARG OPENMS_LIBRARY_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}

RUN apt-get update -y && apt-get install -y python-pip python-dev python-numpy
RUN pip install -U pytest setuptools Cython autowrap pandas
RUN cmake -DCMAKE_PREFIX_PATH="/contrib-build/;/usr/;/usr/local" -DBOOST_USE_STATIC=OFF -DHAS_XSERVER=Off -DPYOPENMS=On ${OPENMS_LIBRARY_DIR}

# make OpenMS library
RUN make pyopenms

# install
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}/pyOpenMS
RUN pip install dist/*.whl
WORKDIR /
ENV PATH="${OPENMS_LIBRARY_BUILD_DIR}/bin/:${PATH}"
LABEL software="OpenMS (library+tools+pyopenms)"
