# syntax=docker/dockerfile:1.3-labs
ARG OPENMS_TAG=Release3.0.0
ARG OPENMS_REPO=https://github.com/OpenMS/OpenMS.git
ARG OPENMS_CONTRIB_REPO=https://github.com/OpenMS/contrib.git
ARG CMAKE_VERSION="3.28.1"
ARG CONTRIB_BUILD_DIR="/contrib-build"
ARG OPENMS_LIBRARY_BUILD_DIR="/openms-build"
ARG OPENMS_LIBRARY_DIR="/OpenMS"
ARG NUM_BUILD_CORES=20
ARG MAKEFLAGS="-j${NUM_BUILD_CORES}"


### Installing OpenMS contributing libraries ###
FROM ubuntu:22.04 as contrib
ARG OPENMS_TAG
ARG OPENMS_CONTRIB_REPO
ARG CMAKE_VERSION
ARG CONTRIB_BUILD_DIR

RUN apt-get -y update \
  && apt-get install -y --no-install-recommends --no-install-suggests \
    g++ \
    autoconf \
    automake \
    patch \
    libtool \
    make \
    git \
    gpg \
    wget \
    ca-certificates \
    curl \
    libboost-date-time1.74-dev \
    libboost-iostreams1.74-dev \
    libboost-regex1.74-dev \
    libboost-math1.74-dev \
    libboost-random1.74-dev \
    libsvm-dev \
    libglpk-dev \
    libzip-dev \
    zlib1g-dev \
    libxerces-c-dev \
    libbz2-dev \
    libomp-dev \
    libhdf5-dev\
    qtbase5-dev \
    libqt5svg5-dev \
    libqt5opengl5-dev \
    libeigen3-dev \
    coinor-libcoinmp-dev \
  && update-ca-certificates

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

WORKDIR /
RUN git clone --branch ${OPENMS_TAG} --single-branch ${OPENMS_CONTRIB_REPO} && rm -rf contrib/.git/
WORKDIR ${CONTRIB_BUILD_DIR}

LABEL base.image="ubuntu:22.04"
LABEL version="3.0"
LABEL software="OpenMS (library)"
LABEL software.version="3.1.0-Ubuntu22.04"
LABEL description="C++ libraries and tools for MS/MS data analysis"
LABEL website="http://www.openms.org/"
LABEL documentation="http://www.openms.org/"
LABEL license="http://www.openms.org/"
LABEL tags="Proteomics"
# to link to repo on github container registry
LABEL org.opencontainers.image.source https://github.com/OpenMS/OpenMS

### Buiding the OpenMS library ###
FROM contrib as library
ARG OPENMS_TAG
ARG OPENMS_REPO
ARG CONTRIB_BUILD_DIR
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

RUN git clone --branch ${OPENMS_TAG} --single-branch ${OPENMS_REPO} ${OPENMS_LIBRARY_DIR}
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN /bin/bash -c "cmake -DCMAKE_BUILD_TYPE='Release' -DCMAKE_PREFIX_PATH='${CONTRIB_BUILD_DIR}/;/usr/;/usr/local' -DBOOST_USE_STATIC=OFF ${OPENMS_LIBRARY_DIR}"
RUN make OpenMS


### Building the OpenMS executables ###
FROM library as executables
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

RUN apt-get update && apt-get install -y --no-install-recommends --no-install-suggests openjdk-17-jdk

WORKDIR ${OPENMS_LIBRARY_DIR}
RUN <<-EOF
    mkdir /thirdparty
    git submodule update --init THIRDPARTY
    cp -r THIRDPARTY/All/* /thirdparty
    cp -r THIRDPARTY/Linux/64bit/* /thirdparty
EOF

ENV PATH="/thirdparty/LuciPHOr2:/thirdparty/MSGFPlus:/thirdparty/Sirius:/thirdparty/ThermoRawFileParser:/thirdparty/Comet:/thirdparty/Fido:/thirdparty/MaRaCluster:/thirdparty/Percolator:/thirdparty/SpectraST:/thirdparty/XTandem:/thirdparty/Sage:${PATH}"
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN make TOPP && rm -rf src doc CMakeFiles
WORKDIR /
ENV PATH="/${OPENMS_LIBRARY_BUILD_DIR}}/bin/:${PATH}"


### Building pyOpenMS ###
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
