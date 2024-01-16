# syntax=docker/dockerfile:1.3-labs
ARG OPENMS_VERSION=3.0.0
ARG OPENMS_REPO=https://github.com/OpenMS/OpenMS.git
ARG OPENMS_LIBRARY_BUILD_DIR="/openms-build"
ARG OPENMS_LIBRARY_DIR="/OpenMS"
ARG CONTRIB_BASE_IMAGE=ghcr.io/openms/contrib:${OPENMS_VERSION}
ARG CONTRIB_BUILD_DIR="/contrib-build"
ARG NUM_BUILD_CORES=20
ARG MAKEFLAGS="-j${NUM_BUILD_CORES}"


### Buiding the OpenMS library ###
FROM ${CONTRIB_BASE_IMAGE} as library
ARG OPENMS_VERSION
ARG OPENMS_REPO
ARG CONTRIB_BUILD_DIR
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG MAKEFLAGS
ARG OPENMS_BRANCH="Release${OPENMS_VERSION}"
ENV MAKEFLAGS="${MAKEFLAGS}"

RUN git clone --branch ${OPENMS_BRANCH} --single-branch ${OPENMS_REPO} ${OPENMS_LIBRARY_DIR}
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN /bin/bash -c "cmake -DCMAKE_BUILD_TYPE='Release' -DCMAKE_PREFIX_PATH='${CONTRIB_BUILD_DIR}/;/usr/;/usr/local' -DBOOST_USE_STATIC=OFF ${OPENMS_LIBRARY_DIR}"
RUN make OpenMS
LABEL software="OpenMS (library)"


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
LABEL software="OpenMS (tools)"


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
LABEL software="OpenMS (library+tools+pyopenms)"
