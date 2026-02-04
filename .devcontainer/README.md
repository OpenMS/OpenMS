# OpenMS Dev Container

This directory contains the VS Code Dev Container configuration for OpenMS development.

## Architecture Support

The dev container now supports both x86_64 (AMD64) and ARM64 (Apple Silicon) architectures.

### Automatic Detection (Recommended)

If you have Docker with BuildKit enabled (Docker 19.03+), the correct architecture should be detected automatically when you open the dev container.

### Manual Configuration for ARM64/Apple Silicon

If automatic detection doesn't work (e.g., you see errors about the wrong architecture or the container fails to start):

1. Open `.devcontainer/devcontainer.json`
2. Modify the `"build"` section to include the arch argument:

```json
"build": {
    "dockerfile": "Dockerfile",
    "args": {
        "TARGETARCH": "arm64"
    }
},
```

3. Rebuild the container (VS Code Command Palette → "Dev Containers: Rebuild Container")

### Enabling BuildKit

If TARGETARCH is not being set automatically, ensure Docker BuildKit is enabled:

```bash
export DOCKER_BUILDKIT=1
```

Or set it permanently in your Docker config (`~/.docker/config.json` or `/etc/docker/daemon.json`):

```json
{
  "features": {
    "buildkit": true
  }
}
```

## Troubleshooting

- **Wrong architecture image**: If you're on ARM64 but getting AMD64 image, manually set `TARGETARCH` as described above
- **Build fails**: Ensure Docker BuildKit is enabled
- **Slow performance on ARM64**: Make sure you're using the ARM64 image, not running AMD64 via emulation

## Images

- AMD64: `ghcr.io/openms/contrib_manylinux_2_34:latest-amd64`
- ARM64: `ghcr.io/openms/contrib_manylinux_2_34:latest-arm64`
