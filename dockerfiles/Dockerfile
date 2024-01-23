# syntax=docker/dockerfile:1.3-labs
ARG OPENMS_REPO=https://github.com/OpenMS/OpenMS.git

# default OpenMS version is "latest" for the current tip of master
# when this image is built from releases it will be "3.1.0" or similar
ARG OPENMS_VERSION=latest

# default OpenMS branch is "master" for the current tip of master
# when this image is built from releases it will be "Release3.1.0" or similar
# NOTE: this is passed to git clone, so tags also work
ARG OPENMS_BRANCH=master
ARG OPENMS_LIBRARY_BUILD_DIR="/openms-build"
ARG OPENMS_LIBRARY_DIR="/OpenMS"
ARG CONTRIB_IMAGE_TAG=latest
ARG CONTRIB_BASE_IMAGE=ghcr.io/openms/contrib:${CONTRIB_IMAGE_TAG}
ARG CONTRIB_BUILD_DIR="/contrib-build"
ARG THIRDPARTY_BRANCH="master"
ARG THIRDPARTY_DIR="/thirdparty"
ARG NUM_BUILD_CORES=8
ARG MAKEFLAGS="-j${NUM_BUILD_CORES}"


### OpenMS base ###
FROM ${CONTRIB_BASE_IMAGE} as base
ARG OPENMS_VERSION

# the main reason we're pulling this in is to make sure that the metadata
# below is consistently formatted for releases
LABEL version="${OPENMS_VERSION}"
LABEL software.version="${OPENMS_VERSION}-Ubuntu22.04"
LABEL org.opencontainers.image.source https://github.com/OpenMS/OpenMS


### OpenMS library ###
FROM base as library
ARG OPENMS_REPO
ARG OPENMS_BRANCH
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG CONTRIB_BUILD_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

RUN git clone --depth=1 --branch=${OPENMS_BRANCH} ${OPENMS_REPO} ${OPENMS_LIBRARY_DIR}
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN cmake \
    -DCMAKE_BUILD_TYPE='Release' \
    -DCMAKE_PREFIX_PATH="${CONTRIB_BUILD_DIR}/;/usr/;/usr/local" \
    -DBOOST_USE_STATIC=OFF \
    ${OPENMS_LIBRARY_DIR}
RUN make OpenMS
LABEL software="OpenMS (library)"


### OpenMS tools build ###
FROM library as tools-build
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG THIRDPARTY_BRANCH
ARG THIRDPARTY_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

RUN apt-get update && apt-get install -y --no-install-recommends --no-install-suggests openjdk-17-jdk

WORKDIR ${OPENMS_LIBRARY_DIR}
RUN <<-EOF
    mkdir /thirdparty
    # note, we are not using the submodule here because the .git folder
    # is not copied to this stage. Instead, we clone the repo directly.
    git clone --branch=${THIRDPARTY_BRANCH} --depth=1 https://github.com/OpenMS/THIRDPARTY.git

    # copying only the relevant binaries and leaving the rest behind
    cp -r THIRDPARTY/All/* ${THIRDPARTY_DIR}
    cp -r THIRDPARTY/Linux/64bit/* ${THIRDPARTY_DIR}

    rm -rf THIRDPARTY
EOF

ENV PATH="${THIRDPARTY_DIR}/LuciPHOr2:${THIRDPARTY_DIR}/MSGFPlus:${THIRDPARTY_DIR}/Sirius:${THIRDPARTY_DIR}/ThermoRawFileParser:${THIRDPARTY_DIR}/Comet:${THIRDPARTY_DIR}/Fido:${THIRDPARTY_DIR}/MaRaCluster:${THIRDPARTY_DIR}/Percolator:${THIRDPARTY_DIR}/SpectraST:${THIRDPARTY_DIR}/XTandem:${THIRDPARTY_DIR}/Sage:${PATH}"
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
# RUN make TOPP && rm -rf src doc CMakeFiles
RUN make TOPP
WORKDIR /
ENV PATH="${OPENMS_LIBRARY_BUILD_DIR}/bin:${PATH}"

### OpenMS tools runtime ###
FROM ubuntu:22.04 AS tools
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG CONTRIB_BUILD_DIR
ARG THIRDPARTY_DIR

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

COPY --from=tools-build ${CONTRIB_BUILD_DIR}/bin ${CONTRIB_BUILD_DIR}/bin
COPY --from=tools-build /thirdparty /thirdparty
COPY --from=tools-build ${OPENMS_LIBRARY_BUILD_DIR}/bin ${OPENMS_LIBRARY_BUILD_DIR}/bin
COPY --from=tools-build ${OPENMS_LIBRARY_BUILD_DIR}/lib ${OPENMS_LIBRARY_BUILD_DIR}/lib
COPY --from=tools-build ${OPENMS_LIBRARY_DIR}/share ${OPENMS_LIBRARY_DIR}/share

ENV PATH="${OPENMS_LIBRARY_BUILD_DIR}/bin:${THIRDPARTY_DIR}/LuciPHOr2:${THIRDPARTY_DIR}/MSGFPlus:${THIRDPARTY_DIR}/Sirius:${THIRDPARTY_DIR}/ThermoRawFileParser:${THIRDPARTY_DIR}/Comet:${THIRDPARTY_DIR}/Fido:${THIRDPARTY_DIR}/MaRaCluster:${THIRDPARTY_DIR}/Percolator:${THIRDPARTY_DIR}/SpectraST:${THIRDPARTY_DIR}/XTandem:${THIRDPARTY_DIR}/Sage:${PATH}"

LABEL software="OpenMS (tools)"


### OpenMS tools runtime test ###
FROM tools AS tools-test
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends --no-install-suggests \
    g++ \
    autoconf \
    automake \
    make \
    gpg \
    wget \
    ca-certificates \
    xvfb \
    curl \
  && update-ca-certificates

RUN <<-EOF
    cmake_ubuntu_version=$(lsb_release -cs)
    if ! wget -q --method=HEAD "https://apt.kitware.com/ubuntu/dists/$cmake_ubuntu_version/Release"; then
      bash -c "$(wget -O - https://apt.kitware.com/kitware-archive.sh)"
    else
      wget -qO - https://apt.kitware.com/kitware-archive.sh | bash -s -- --release $cmake_ubuntu_version
    fi
    apt-get -y update
    apt-get install -y cmake
EOF

COPY --from=tools-build ${OPENMS_LIBRARY_BUILD_DIR} ${OPENMS_LIBRARY_BUILD_DIR}
COPY --from=tools-build ${OPENMS_LIBRARY_DIR}/src/tests ${OPENMS_LIBRARY_DIR}/src/tests
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN <<-EOF
    # we need Xvfb to run a small subset of tests
    # mostly TOPP_INIUpdater
    Xvfb :0 &
    DISPLAY=:0 ctest -R TOPP_
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
