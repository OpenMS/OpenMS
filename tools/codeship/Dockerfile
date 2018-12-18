FROM openms/contrib:qt5-v2.3.1

RUN mkdir /OpenMS
WORKDIR /OpenMS
COPY . .

WORKDIR /
RUN mkdir openms-build
WORKDIR /openms-build

# Load the Qt5 environment and configure
RUN QT_ENV=$(find /opt -name 'qt*-env.sh') && /bin/bash -c "source ${QT_ENV} && cmake -DCMAKE_PREFIX_PATH='/contrib-build/;/usr/;/usr/local' -DBOOST_USE_STATIC=OFF ../OpenMS"
RUN make OpenMS -j2

