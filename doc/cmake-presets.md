# CMake Presets for OpenMS

This document describes the CMake presets available for building OpenMS. CMake presets provide a standardized way to configure common build scenarios.

## What are CMake Presets?

CMake presets (introduced in CMake 3.19, fully featured in 3.21+) allow you to specify common CMake configuration options in a JSON file. This simplifies the build process and ensures consistency across different build environments.

For more information, see the [official CMake presets documentation](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html).

## Available Presets

### Basic Build Types

- **`release`**: Standard release build with optimizations
- **`debug`**: Debug build with symbols
- **`relwithdebinfo`**: Release build with debug information

### CI-Specific Presets

These presets are optimized for continuous integration builds:

- **`ci-linux`**: Linux CI build with GCC
- **`ci-linux-clang`**: Linux CI build with Clang
- **`ci-macos`**: macOS CI build
- **`ci-windows`**: Windows CI build with Visual Studio 2022

CI presets disable git tracking, update checks, and tutorials by default.

### pyOpenMS Presets

For building the Python bindings:

- **`pyopenms-linux`**: pyOpenMS build for Linux
- **`pyopenms-macos`**: pyOpenMS build for macOS
- **`pyopenms-windows`**: pyOpenMS build for Windows

### Package Building Presets

For creating distribution packages:

- **`package-deb`**: Debian/Ubuntu package
- **`package-rpm`**: RPM package (Fedora, RHEL, etc.)
- **`package-dmg`**: macOS DMG image
- **`package-pkg`**: macOS PKG installer
- **`package-nsis`**: Windows NSIS installer
- **`package-zip`**: ZIP archive (cross-platform)

### KNIME Plugin Presets

For building KNIME packages:

- **`knime-linux`**: KNIME package for Linux
- **`knime-macos`**: KNIME package for macOS
- **`knime-windows`**: KNIME package for Windows

### Library Dependency Presets

- **`system-libs`**: Use system-installed libraries
- **`contrib-libs`**: Use contrib libraries (requires setting OPENMS_CONTRIB_LIBS)

## Using Presets

### List Available Presets

```bash
cmake --list-presets
```

### Configure with a Preset

```bash
cmake --preset=<preset-name>
```

For example:
```bash
cmake --preset=release
```

### Build with a Preset

```bash
cmake --build --preset=<preset-name>
```

For example:
```bash
cmake --build --preset=pyopenms-linux
```

### Test with a Preset

```bash
ctest --preset=<preset-name>
```

For example:
```bash
ctest --preset=ci-linux
```

## Creating Custom User Presets

You can create your own presets by creating a `CMakeUserPresets.json` file in the repository root. This file is git-ignored and won't be committed.

See `CMakeUserPresets.json.example` for examples of how to create custom presets.

Example workflow for custom preset:

1. Copy the example file:
   ```bash
   cp CMakeUserPresets.json.example CMakeUserPresets.json
   ```

2. Edit `CMakeUserPresets.json` to match your system paths

3. Use your custom preset:
   ```bash
   cmake --preset=user-example
   ```

## Common Customizations

### Setting Contrib Library Path

Most presets assume you have a contrib build. Set the path using environment variables or in your user presets:

```json
{
  "configurePresets": [
    {
      "name": "my-build",
      "inherits": "release",
      "cacheVariables": {
        "OPENMS_CONTRIB_LIBS": "/path/to/your/contrib"
      }
    }
  ]
}
```

### Setting Qt Path

If Qt is not in your system path:

```json
{
  "configurePresets": [
    {
      "name": "my-build",
      "inherits": "release",
      "cacheVariables": {
        "CMAKE_PREFIX_PATH": "/path/to/Qt/6.x.x/gcc_64"
      }
    }
  ]
}
```

### Specifying Python for pyOpenMS

```json
{
  "configurePresets": [
    {
      "name": "my-pyopenms",
      "inherits": "pyopenms-linux",
      "cacheVariables": {
        "Python_EXECUTABLE": "/path/to/python3"
      }
    }
  ]
}
```

## Using with IDEs

Many modern IDEs support CMake presets:

- **Visual Studio 2022**: Opens `CMakePresets.json` automatically
- **Visual Studio Code**: Use the CMake Tools extension
- **CLion**: Supports presets starting from version 2021.3

## Environment Variables

You can use environment variables in preset configurations:

```json
{
  "configurePresets": [
    {
      "name": "my-preset",
      "cacheVariables": {
        "OPENMS_CONTRIB_LIBS": "$env{CONTRIB_PATH}"
      }
    }
  ]
}
```

Then set the environment variable before running cmake:
```bash
export CONTRIB_PATH=/path/to/contrib
cmake --preset=my-preset
```

## Preset Inheritance

Presets use inheritance to avoid duplication. For example:

- `ci-linux` inherits from `ci-base` which inherits from `base`
- You can create your own presets that inherit from existing ones

## Migrating from Shell Scripts

If you're currently using shell scripts for configuration (like in CI), presets can replace many of those configurations. For example, instead of:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DPYOPENMS=ON -DWITH_GUI=OFF ..
```

You can use:
```bash
cmake --preset=pyopenms-linux
```

## Benefits of Using Presets

1. **Consistency**: Same configuration across different developers and CI systems
2. **Simplicity**: Single command instead of long cmake invocations
3. **Documentation**: Presets serve as documentation of common configurations
4. **IDE Support**: Better integration with modern IDEs
5. **Inheritance**: Share common settings across multiple configurations
6. **Version Control**: Presets are versioned with the code

## Troubleshooting

### Preset Not Found

Make sure you're using CMake 3.21 or later:
```bash
cmake --version
```

### Preset Conditions Not Met

Some presets have conditions (e.g., only available on certain platforms). Check the preset description and conditions in `CMakePresets.json`.

### Conflicting Cache Variables

If you've previously configured the build directory with different settings, try:
```bash
rm -rf bld
cmake --preset=<preset-name>
```

## Additional Resources

- [CMake Presets Official Documentation](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)
- [OpenMS Build Documentation](https://openms.readthedocs.io)
