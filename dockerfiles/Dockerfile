# syntax=docker/dockerfile:1.3-labs
# default OpenMS version is "latest" for the current tip of master
# when this image is built from releases it will be "3.1.0" or similar
ARG OPENMS_VERSION_TAG=latest
ARG SOURCE_DIR="/tmp/OpenMS"
ARG BUILD_DIR="${SOURCE_DIR}/bld"
ARG INSTALL_DIR="/opt/OpenMS"
ARG THIRDPARTY_DIR="/opt/OpenMS/thirdparty"
ARG NUM_BUILD_CORES=8
ARG MAKEFLAGS="-j${NUM_BUILD_CORES}"


### OpenMS base ###
FROM ubuntu:24.04 AS base
ARG OPENMS_VERSION_TAG
ARG CMAKE_VERSION=3.28.1
ARG TARGETARCH

RUN apt-get update \
    && apt-get install -y software-properties-common \
    && add-apt-repository universe \
    && apt-get update \
    && apt-get install -y --no-install-recommends --no-install-suggests \
    g++ \
    make \
    git \
    ca-certificates \
    build-essential \
    # OpenMS build dependencies
    libsvm-dev \
    libglpk-dev \
    libzip-dev \
    zlib1g-dev \
    libzip4 \
    zlib1g \
    libxerces-c-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libomp-dev \
    libhdf5-dev\
    coinor-libcoinmp-dev \
    libeigen3-dev \
    libboost-date-time-dev \
    libboost-iostreams-dev \
    libboost-regex-dev \
    libboost-math-dev \
    libboost-random-dev \
    qt6-base-dev \
    libglx-dev \
    libgl1-mesa-dev \
    libqt6opengl6-dev \
    libqt6svg6-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && apt-get update && apt-get install -y --no-install-recommends wget \
  && wget https://repo1.maven.org/maven2/org/apache/arrow/ubuntu/apache-arrow-apt-source-latest-noble.deb \
  && apt-get install -y --no-install-recommends ./apache-arrow-apt-source-latest-noble.deb \
  # Pin Arrow/Parquet to v23 so build-time dev libs match the runtime SONAME (libarrow2300/libparquet2300)
  && printf 'Package: libarrow* libparquet*\nPin: version 23.*\nPin-Priority: 999\n' > /etc/apt/preferences.d/arrow \
  && apt-get update && apt-get install -y --no-install-recommends --no-install-suggests libparquet-dev \
  && rm apache-arrow-apt-source-latest-noble.deb

# installing cmake
WORKDIR /tmp
RUN <<-EOF
    set -eux
    # Map Docker's TARGETARCH to CMake's architecture naming
    case "${TARGETARCH}" in
        amd64) CMAKE_ARCH="x86_64" ;;
        arm64) CMAKE_ARCH="aarch64" ;;
        *) echo "Unsupported architecture: ${TARGETARCH}" && exit 1 ;;
    esac
    wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-${CMAKE_ARCH}.sh -O cmake.sh
    mkdir -p /opt/cmake
    sh cmake.sh --skip-license --prefix=/opt/cmake
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
    ln -s /opt/cmake/bin/ctest /usr/local/bin/ctest
    rm -rf /tmp/*
EOF

# The main reason we're pulling this in is to make sure that the metadata
# below is consistently formatted for releases irrespective of the base image
LABEL version="${OPENMS_VERSION_TAG}"
LABEL software.version="${OPENMS_VERSION_TAG}-Ubuntu24.04"
LABEL org.opencontainers.image.source https://github.com/OpenMS/OpenMS


### OpenMS build ###
FROM base AS build
ARG SOURCE_DIR
ARG BUILD_DIR
ARG INSTALL_DIR
ARG THIRDPARTY_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

COPY . ${SOURCE_DIR}

WORKDIR ${SOURCE_DIR}
RUN mkdir -p ${THIRDPARTY_DIR}
RUN <<-EOF
    set -eux

    # Map Docker's TARGETARCH to THIRDPARTY directory naming
    case "${TARGETARCH}" in
        amd64) ARCH_DIR="x86_64" ;;
        arm64) ARCH_DIR="aarch64" ;;
        *) echo "Unsupported architecture: ${TARGETARCH}" && exit 1 ;;
    esac

    ALL_DIR=THIRDPARTY/All
    LINUX_DIR=THIRDPARTY/Linux/${ARCH_DIR}
    if [ ! -d "$ALL_DIR" ] || [ ! -d "$LINUX_DIR" ]; then
      git submodule update --depth=1 --init THIRDPARTY
    fi

    # copying only the binaries that are relevant to Linux
    # Copy All directory (platform-independent files)
    if [ -d "$ALL_DIR" ]; then
      cp -vr "$ALL_DIR"/* ${THIRDPARTY_DIR}
    fi
    
    # Copy architecture-specific binaries if they exist
    if [ -d "$LINUX_DIR" ]; then
      cp -vr "$LINUX_DIR"/* ${THIRDPARTY_DIR}
    else
      echo "Warning: No architecture-specific binaries found at $LINUX_DIR"
      echo "Building without third-party tools for ${ARCH_DIR}"
    fi
    
    rm -rf THIRDPARTY
EOF

# adding thirdparty binaries to the PATH so they will be picked up by the build system
ENV PATH="${THIRDPARTY_DIR}/LuciPHOr2:${THIRDPARTY_DIR}/MSGFPlus:${THIRDPARTY_DIR}/ThermoRawFileParser:${THIRDPARTY_DIR}/Comet:${THIRDPARTY_DIR}/Fido:${THIRDPARTY_DIR}/MaRaCluster:${THIRDPARTY_DIR}/Percolator:${THIRDPARTY_DIR}/SpectraST:${THIRDPARTY_DIR}/Sage:${PATH}"
RUN <<-EOF
    set -eux
    # Map Docker's TARGETARCH to Debian multiarch tuple
    case "${TARGETARCH}" in
        amd64) DEB_ARCH="x86_64-linux-gnu" ;;
        arm64) DEB_ARCH="aarch64-linux-gnu" ;;
        *) echo "Unsupported architecture: ${TARGETARCH}" && exit 1 ;;
    esac
    
    cmake \
        -DQt6_DIR="/usr/lib/${DEB_ARCH}/cmake/Qt6" \
        -DCMAKE_BUILD_TYPE='Release' \
        -DCMAKE_PREFIX_PATH="/usr/;/usr/local" \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DBOOST_USE_STATIC=OFF \
        -S ${SOURCE_DIR} \
        -B ${BUILD_DIR}
EOF

WORKDIR ${BUILD_DIR}
RUN make all
RUN make install/strip


### OpenMS library ###
FROM ubuntu:24.04 AS library
ARG INSTALL_DIR

RUN apt-get update \
  && apt-get install -y --no-install-recommends --no-install-suggests \
    ca-certificates \
    wget \
    qt6-base-dev \
    libqt6opengl6 \
    libsvm3 \
    libzip4 \
    zlib1g \
    libbz2-1.0 \
    libcurl4 \
    libgomp1 \
    libqt6svg6 \
    libxerces-c3.2 \
    coinor-libcoinmp1v5 \
    libboost-date-time1.83.0 \
    libboost-iostreams1.83.0 \
    libboost-regex1.83.0 \
    libboost-math1.83.0 \
    libboost-random1.83.0 \
  && wget https://repo1.maven.org/maven2/org/apache/arrow/ubuntu/apache-arrow-apt-source-latest-noble.deb \
  && apt-get install -y --no-install-recommends ./apache-arrow-apt-source-latest-noble.deb \
  && apt-get update && apt-get install -y --no-install-recommends --no-install-suggests \
    libarrow2300 \
    libparquet2300 \
  && rm apache-arrow-apt-source-latest-noble.deb \
  && apt-get purge -y wget && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/*

COPY --from=build ${INSTALL_DIR}/lib ${INSTALL_DIR}/lib
COPY --from=build ${INSTALL_DIR}/include ${INSTALL_DIR}/include
COPY --from=build ${INSTALL_DIR}/share ${INSTALL_DIR}/share

LABEL software="OpenMS (library)"


### OpenMS tools runtime ###
FROM library AS tools
ARG INSTALL_DIR

COPY --from=build ${INSTALL_DIR}/bin ${INSTALL_DIR}/bin
ENV PATH="${INSTALL_DIR}/bin:${PATH}"

# Ensure a writable home directory for config files (e.g., UpdateCheck .ver files).
# Without this, $HOME may resolve to / in CI environments, causing tools to crash
# when trying to create $HOME/.config/OpenMS.
# This ENV propagates to all child stages (tools-thirdparty, test).
RUN mkdir -p /app && chmod -R 755 /app
ENV HOME=/app

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

ENV PATH="${THIRDPARTY_DIR}/LuciPHOr2:${THIRDPARTY_DIR}/MSGFPlus:${THIRDPARTY_DIR}/ThermoRawFileParser:${THIRDPARTY_DIR}/Comet:${THIRDPARTY_DIR}/Fido:${THIRDPARTY_DIR}/MaRaCluster:${THIRDPARTY_DIR}/Percolator:${THIRDPARTY_DIR}/SpectraST:${THIRDPARTY_DIR}/Sage:${PATH}"

LABEL software="OpenMS (tools + thirdparty)"


### OpenMS runtime test ###
# this stage is based on the runtime stage, so it's a test to ensure
# that the runtime dependencies satisfy the requirements of the tools,
# binaries and the OpenMS library
FROM tools-thirdparty AS test
ARG SOURCE_DIR
ARG BUILD_DIR
ARG NUM_BUILD_CORES=8
ARG CMAKE_VERSION=3.28.1
ARG TARGETARCH

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends --no-install-suggests \
    # we need Xvfb to run a small subset of tests (eg. TOPP_INIUpdater)
    xvfb \
    # needed for TSGDialog_test and TOPPView_test
    libqt6test6 \
    wget \
  && rm -rf /var/lib/apt/lists/*

# installing cmake
WORKDIR /tmp
RUN <<-EOF
    set -eux
    # Map Docker's TARGETARCH to CMake's architecture naming
    case "${TARGETARCH}" in
        amd64) CMAKE_ARCH="x86_64" ;;
        arm64) CMAKE_ARCH="aarch64" ;;
        *) echo "Unsupported architecture: ${TARGETARCH}" && exit 1 ;;
    esac
    wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-${CMAKE_ARCH}.sh -O cmake.sh
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
