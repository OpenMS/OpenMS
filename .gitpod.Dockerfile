FROM openms/contrib:latest
                    
#USER gitpod
# Avoid user input
ARG DEBIAN_FRONTEND=noninteractive

# Install tools for VSCode, Intellisense, JRE for Thirdparties, etc. using apt-get
RUN apt-get -q update && apt-get install -yq gdb unzip wget php openjdk-11-jre clang-tidy && rm -rf /var/lib/apt/lists/*
RUN wget https://github.com/clangd/clangd/releases/download/12.0.1/clangd-linux-12.0.1.zip && unzip clangd-linux-12.0.1.zip -d /usr/ && rm clangd-linux-12.0.1.zip 
#
# More information: https://www.gitpod.io/docs/42_config_docker/
