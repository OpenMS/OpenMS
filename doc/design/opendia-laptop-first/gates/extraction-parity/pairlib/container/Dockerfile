FROM ubuntu:24.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libgomp1 libicu74 zstd time ca-certificates python3 && rm -rf /var/lib/apt/lists/*
ADD payload.tar /
ENV DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
ENV PATH=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet:/usr/local/bin:/usr/bin:/bin
ENV OPENMS_DATA_PATH=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/share/OpenMS
