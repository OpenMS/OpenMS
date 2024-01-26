# syntax=docker/dockerfile:1.3-labs
ARG OPENMS_REPO=https://github.com/OpenMS/OpenMS.git

# default OpenMS version is "latest" for the current tip of master
# when this image is built from releases it will be "3.1.0" or similar
ARG OPENMS_VERSION=latest

# default OpenMS branch is "master" for the current tip of master
# when this image is built from releases it will be "Release3.1.0" or similar
# NOTE: this is passed to git clone, so tags also work
ARG OPENMS_BRANCH=master

# this is the directory where OPENMS_REPO gets cloned to
ARG OPENMS_SOURCE_DIR="/tmp/OpenMS"

# this is the directory where the OpenMS library and tools are built
ARG OPENMS_LIBRARY_BUILD_DIR="/opt/OpenMS"
ARG OPENMS_LIBRARY_DIR="/opt/OpenMS/lib"
ARG CONTRIB_IMAGE_TAG=latest
ARG CONTRIB_BASE_IMAGE=ghcr.io/openms/contrib:${CONTRIB_IMAGE_TAG}
ARG CONTRIB_BUILD_DIR="/contrib-build"
ARG THIRDPARTY_BRANCH="master"
ARG THIRDPARTY_DIR="/opt/OpenMS/thirdparty"
ARG NUM_BUILD_CORES=8
ARG MAKEFLAGS="-j${NUM_BUILD_CORES}"
ARG CMAKE_VERSION="3.28.1"


### OpenMS base ###
FROM ${CONTRIB_BASE_IMAGE} as base
ARG OPENMS_VERSION

# the main reason we're pulling this in is to make sure that the metadata
# below is consistently formatted for releases
LABEL version="${OPENMS_VERSION}"
LABEL software.version="${OPENMS_VERSION}-Ubuntu22.04"
LABEL org.opencontainers.image.source https://github.com/OpenMS/OpenMS


### OpenMS library build ###
FROM base as library-build
ARG OPENMS_REPO
ARG OPENMS_SOURCE_DIR
ARG OPENMS_BRANCH
ARG OPENMS_LIBRARY_BUILD_DIR
ARG CONTRIB_BUILD_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

RUN git clone --depth=1 --branch=${OPENMS_BRANCH} ${OPENMS_REPO} ${OPENMS_SOURCE_DIR}
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
RUN <<-EOF
    git submodule update --init THIRDPARTY --depth=1

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
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG CONTRIB_BUILD_DIR
ARG THIRDPARTY_DIR

RUN apt-get update \
  && apt-get install -y --no-install-recommends --no-install-suggests \
    ## thirdparty dependencies
    # LuciPHOr2 and MSGFPlus:
    openjdk-17-jre-headless \
  && rm -rf /var/lib/apt/lists/*

COPY --from=tools-build ${OPENMS_LIBRARY_BUILD_DIR}/bin ${OPENMS_LIBRARY_BUILD_DIR}/bin

ENV PATH="${OPENMS_LIBRARY_BUILD_DIR}/bin:${THIRDPARTY_DIR}/LuciPHOr2:${THIRDPARTY_DIR}/MSGFPlus:${THIRDPARTY_DIR}/ThermoRawFileParser:${THIRDPARTY_DIR}/Comet:${THIRDPARTY_DIR}/Fido:${THIRDPARTY_DIR}/MaRaCluster:${THIRDPARTY_DIR}/Percolator:${THIRDPARTY_DIR}/SpectraST:${THIRDPARTY_DIR}/XTandem:${THIRDPARTY_DIR}/Sage:${PATH}"
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}

LABEL software="OpenMS (tools)"


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

# installing cmake
WORKDIR /tmp
ADD https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.sh cmake.sh
RUN <<-EOF
    mkdir -p /opt/cmake
    sh cmake.sh --skip-license --prefix=/opt/cmake
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
    ln -s /opt/cmake/bin/ctest /usr/local/bin/ctest
    rm -rf /tmp/*
EOF

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends --no-install-suggests \
    # we need Xvfb to run a small subset of tests (eg. TOPP_INIUpdater)
    xvfb \
  && rm -rf /var/lib/apt/lists/*

COPY --from=test-build ${OPENMS_SOURCE_DIR} ${OPENMS_SOURCE_DIR}
COPY --from=test-build ${OPENMS_LIBRARY_BUILD_DIR} ${OPENMS_LIBRARY_BUILD_DIR}
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN <<-EOF
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
