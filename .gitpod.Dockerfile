FROM ghcr.io/openms/contrib
                    
#USER gitpod
# Avoid user input
ARG DEBIAN_FRONTEND=noninteractive
ARG CLANGDVER=14.0.0

# Install tools for VSCode, Intellisense, JRE for Thirdparties, etc. using apt-get
RUN apt-get -q update && apt-get install -yq gdb zip tar unzip curl wget php openjdk-11-jre mono-complete python3-pip && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/clangd/clangd/releases/download/$CLANGDVER/clangd-linux-$CLANGDVER.zip && unzip clangd-linux-$CLANGDVER.zip -d /opt/ && mv /opt/clangd_$CLANGDVER /opt/clangd/ && rm clangd-linux-$CLANGDVER.zip 
#
# More information: https://www.gitpod.io/docs/42_config_docker/
