FROM openms/contrib:latest
                    
#USER gitpod

# Install custom tools, runtime, etc. using apt-get
RUN sudo apt-get -q update && sudo apt-get install -yq openjdk-8-jre && sudo rm -rf /var/lib/apt/lists/*
#
# More information: https://www.gitpod.io/docs/42_config_docker/
