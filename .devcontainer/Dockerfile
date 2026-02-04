# Select the correct base image based on architecture
# TARGETARCH should be set by Docker BuildKit, but we provide a default
# and fallback logic for environments where BuildKit is not enabled

ARG TARGETARCH

# For environments without BuildKit, detect architecture at runtime
FROM alpine:latest AS archdetect
RUN if [ -z "${TARGETARCH}" ]; then \
      ARCH=$(uname -m); \
      case "$ARCH" in \
        x86_64) echo "amd64" > /tmp/arch.txt ;; \
        aarch64|arm64) echo "arm64" > /tmp/arch.txt ;; \
        *) echo "amd64" > /tmp/arch.txt ;; \
      esac; \
    else \
      echo "${TARGETARCH}" > /tmp/arch.txt; \
    fi && cat /tmp/arch.txt

# Unfortunately, we cannot use the detected arch in FROM dynamically,
# so we need TARGETARCH to be set. Default to amd64 if not provided.
FROM ghcr.io/openms/contrib_manylinux_2_34:latest-${TARGETARCH:-amd64}

# Copy detection result for debugging
COPY --from=archdetect /tmp/arch.txt /tmp/detected_arch.txt
RUN echo "Target architecture: ${TARGETARCH:-amd64}" && \
    echo "Detected architecture: $(cat /tmp/detected_arch.txt 2>/dev/null || echo 'unknown')"

