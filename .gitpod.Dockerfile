FROM openms/contrib:latest
                    
#USER gitpod

# Install custom tools, runtime, etc. using apt-get
RUN apt-get -q update && apt-get install -yq openjdk-8-jre && rm -rf /var/lib/apt/lists/*
#
# More information: https://www.gitpod.io/docs/42_config_docker/
