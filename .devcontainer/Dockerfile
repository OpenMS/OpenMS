# Select the correct base image based on architecture
ARG TARGETARCH
FROM ghcr.io/openms/contrib_manylinux_2_34:latest-${TARGETARCH}

# TARGETARCH is automatically set by Docker BuildKit to the target architecture
# It will be "amd64" on x86_64 systems and "arm64" on ARM64 systems
