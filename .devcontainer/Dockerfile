# TARGETARCH is set automatically by Docker BuildKit (default since Docker 23.0).
# Falls back to amd64 when BuildKit is not available.
# ARM64/Apple Silicon users: if auto-detection fails, add to devcontainer.json:
#   "build": { "args": { "TARGETARCH": "arm64" } }
ARG TARGETARCH=amd64
FROM ghcr.io/openms/contrib_manylinux_2_34:latest-${TARGETARCH}
