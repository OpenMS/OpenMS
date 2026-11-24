FROM ubuntu:24.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libgomp1 libicu74 ca-certificates time && rm -rf /var/lib/apt/lists/*
ADD payload.tar /
ENV DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
ENV PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin:/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet:/usr/local/bin:/usr/bin:/bin
ENV OPENMS_DATA_PATH=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/share/OpenMS
ENV LD_LIBRARY_PATH=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib
