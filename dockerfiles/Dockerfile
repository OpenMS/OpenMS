# syntax=docker/dockerfile:1.3-labs
ARG OPENMS_VERSION=3.1.0
ARG OPENMS_LIBRARY_BUILD_DIR="/openms-build"
ARG OPENMS_LIBRARY_DIR="/OpenMS"
ARG CONTRIB_TAG=latest
ARG CONTRIB_BASE_IMAGE=ghcr.io/openms/contrib:${CONTRIB_TAG}
ARG CONTRIB_BUILD_DIR="/contrib-build"
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
ARG CONTRIB_BUILD_DIR
ARG OPENMS_LIBRARY_DIR
ARG OPENMS_LIBRARY_BUILD_DIR
ARG MAKEFLAGS
ENV MAKEFLAGS="${MAKEFLAGS}"

# NOTE: this assumes that the context is root of the OpenMS repo
COPY . ${OPENMS_LIBRARY_DIR}
WORKDIR ${OPENMS_LIBRARY_BUILD_DIR}
RUN /bin/bash -c "cmake -DCMAKE_BUILD_TYPE='Release' -DCMAKE_PREFIX_PATH='${CONTRIB_BUILD_DIR}/;/usr/;/usr/local' -DBOOST_USE_STATIC=OFF ${OPENMS_LIBRARY_DIR}"
RUN make OpenMS
LABEL software="OpenMS (library)"


### OpenMS tools ###
FROM library as tools
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
