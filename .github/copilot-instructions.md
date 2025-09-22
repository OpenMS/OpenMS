# GitHub Copilot Instructions for OpenMS

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Project Overview

OpenMS is an open-source software C++ library for LC-MS data management and analyses. It offers an infrastructure for rapid development of mass spectrometry related software. OpenMS is free software available under the three-clause BSD license and runs under Windows, macOS, and Linux.

### Key Features

- **Core C++ library** with modern C++20 standards
- **Python bindings** (pyOpenMS) for rapid algorithm development
- **150+ command-line tools** (TOPP Tools) for MS data processing
- **Visualization tools** (TOPPView) for 1D, 2D, and 3D data
- **Support for major file formats** (mzML, mzXML, mzIdentML, pepXML, mzTab)
- **Comprehensive quantitation support** (label-free, SILAC, iTRAQ, TMT, SRM, SWATH)
- **Integration with workflow systems** (KNIME, Galaxy, nextflow)

## Working Effectively

Bootstrap, build, and test the repository:

### Dependencies Installation (Ubuntu/Linux - takes 3-5 minutes)
```bash
# Update package repositories
sudo add-apt-repository universe -y && sudo apt update

# Install required build dependencies
sudo apt-get -qq install -y \
  build-essential cmake autoconf patch libtool git automake ninja-build xvfb ccache \
  qt6-base-dev libqt6svg6-dev libqt6opengl6-dev libqt6openglwidgets6 libgl-dev \
  libeigen3-dev libboost-random-dev libboost-regex-dev libboost-iostreams-dev \
  libboost-date-time-dev libboost-math-dev libxerces-c-dev zlib1g-dev libsvm-dev \
  libbz2-dev coinor-libcoinmp-dev libhdf5-dev

# Optional dependencies for documentation
sudo apt-get -qq install -y doxygen ghostscript graphviz
```

### Contrib Dependencies Build (NEVER CANCEL - takes 10-15 minutes)
```bash
# Initialize contrib submodule
git submodule update --init contrib

# Build contrib libraries - NEVER CANCEL: Set timeout to 20+ minutes
mkdir -p contrib-build
cd contrib-build
cmake -DBUILD_TYPE=ALL -DNUMBER_OF_JOBS=4 ../contrib
make -j 4
cd ..
```

### Main OpenMS Build (NEVER CANCEL - takes 45-90 minutes)  
```bash
# Build OpenMS library and tools - NEVER CANCEL: Set timeout to 120+ minutes
mkdir -p openms-build
cd openms-build
CONTRIB_PATH=`pwd`/../contrib-build
cmake -DOPENMS_CONTRIB_LIBS="$CONTRIB_PATH" -DBOOST_USE_STATIC=On ..
make -j 4
cd ..
```

### Testing (NEVER CANCEL - takes 20-45 minutes)
```bash
# Run unit tests - NEVER CANCEL: Set timeout to 60+ minutes
cd openms-build
ctest -j 4

# Run specific test categories
ctest -R "class_tests" -j 4      # Unit tests only
ctest -R "topp" -j 4             # TOPP tool tests only
ctest -R "pyopenms" -j 4         # Python binding tests only
```

### Quick Build Script Alternative
Use the provided quick build script:
```bash
# Quick build with default settings - NEVER CANCEL: Set timeout to 120+ minutes
./tools/quickbuild.sh 4  # Use 4 parallel jobs
```

## Validation

### Build Validation
Always verify the build completed successfully:
```bash
# Check that key components were built
ls openms-build/bin/       # Should contain TOPP tools
ls openms-build/lib/       # Should contain OpenMS libraries
ls contrib-build/lib/      # Should contain contrib libraries
```

### TOPP Tools Validation
Test that TOPP tools work correctly:
```bash
cd openms-build
# Test basic TOPP tool functionality
./bin/FileInfo --help
./bin/FeatureFinderCentroided --help
./bin/PeptideIndexer --help
```

### Manual Testing Scenarios
After making changes, always run through these validation scenarios:

1. **Basic Library Functionality**: Verify core data structures compile and link correctly
2. **TOPP Tool Integration**: Test that modified algorithms work in command-line tools  
3. **Python Bindings**: If touching core classes, verify pyOpenMS bindings still work
4. **Cross-platform Compatibility**: Consider impact on Windows and macOS builds

## Build System Architecture

### Key Build Targets
```bash
# Available CMake targets (from cmake/messages.cmake):
make OpenMS              # Core library only
make TOPP               # All TOPP command-line tools  
make GUI                # GUI tools (TOPPView, TOPPAS)
make test               # Run all tests
make doc                # Generate documentation
make pyopenms           # Python bindings
```

### Important Build Variables
- `OPENMS_CONTRIB_LIBS`: Path to built contrib libraries (required)
- `BOOST_USE_STATIC`: Use static Boost libraries (recommended: On)
- `CMAKE_BUILD_TYPE`: Release, Debug, RelWithDebInfo (default: Release)
- `ENABLE_UPDATE_CHECK`: Update notifications (default: On)

## Project Structure

```
OpenMS/
├── cmake/                   # CMake build system files
├── contrib/                 # Third-party dependencies (submodule)
├── doc/                     # Documentation
├── src/
│   ├── openms/              # Core library
│   │   ├── include/OpenMS/  # Header files
│   │   └── source/          # Implementation files
│   ├── openms_gui/          # GUI components
│   ├── pyOpenMS/            # Python bindings
│   │   ├── pxds/            # Declarations for autowrap
│   │   └── addons/          # Manual wrapper code
│   ├── tests/               # Test suites
│   │   ├── class_tests/     # Unit tests
│   │   └── topp/            # TOPP tool tests
│   └── topp/                # TOPP tools
├── tools/                   # Development utilities
│   ├── ci/                  # CI scripts
│   └── quickbuild.sh        # Quick build script
└── .github/workflows/       # GitHub Actions CI
```

## Testing Infrastructure

### Test Organization
- **Class tests**: Unit tests in `src/tests/class_tests/openms/`
- **TOPP tests**: Integration tests in `src/tests/topp/`  
- **Python tests**: pyOpenMS tests in `src/pyOpenMS/tests/`
- Follow naming convention: `ClassNameTest.cpp` for `ClassName`

### Running Specific Tests
```bash
cd openms-build
# Run tests for specific class
ctest -R "ClassName_test"

# Run tests with verbose output
ctest -V -R "test_pattern"

# Run tests matching pattern
ctest -R "TOPP.*Tool.*test"
```

## Development Guidelines

### Code Style and Conventions
- Follow existing C++ coding conventions in the codebase
- Use established naming patterns for classes, methods, and variables
- Add Doxygen comments for new public methods and classes
- Maintain consistency with surrounding code style

### Common Development Tasks
- Always check existing tests pass before making changes
- Write unit tests for new functionality
- Use OpenMS data structures (MSExperiment, FeatureMap, PeptideIdentification)
- Follow established error handling patterns
- Consider memory usage when processing large datasets

### Pre-commit Validation
Always run these commands before committing changes:
```bash
# Build with your changes
make -j 4

# Run relevant tests  
ctest -R "affected_component"

# Check code formatting (if available)
# clang-format validation happens in CI
```

## Common Issues and Solutions

### Build Issues
- **Out of memory**: Reduce parallel jobs (`make -j 2` instead of `make -j 4`)
- **Missing contrib**: Ensure `git submodule update --init contrib` was run
- **Qt6 issues**: Verify qt6-base-dev and related packages are installed

### Test Issues
- **Test failures**: Check if tests were already failing before your changes
- **Timeout issues**: Some tests require large timeouts due to computational complexity
- **Missing test data**: Some tests require data files in specific locations

## CI Integration

The project uses GitHub Actions for CI (`.github/workflows/`):
- **openms_ci_matrix_full.yml**: Main CI pipeline for multiple platforms
- **test_pyopenms.yml**: Python binding specific tests
- Tests run on Ubuntu, macOS, and Windows

## Time Expectations

**CRITICAL TIMING INFORMATION - NEVER CANCEL THESE OPERATIONS:**

- **Dependencies installation**: 3-5 minutes
- **Contrib build**: 10-15 minutes  
- **Main OpenMS build**: 45-90 minutes
- **Full test suite**: 20-45 minutes
- **Documentation build**: 10-20 minutes

Always set timeouts with significant buffers:
- Contrib build: 30+ minute timeout
- Main build: 120+ minute timeout  
- Test execution: 60+ minute timeout

## Memory and Resource Requirements

- **Minimum RAM**: 8GB (16GB recommended for parallel builds)
- **Disk space**: 10GB+ for full build with contrib
- **CPU cores**: Build scales with available cores (tested up to 4 cores)

## Assistance Focus

When providing code suggestions:
- Analyze existing code patterns before suggesting changes
- Respect the established architecture and build system
- Consider cross-platform compatibility (Windows, Linux, macOS)
- Use existing OpenMS data structures and algorithms
- Follow established error handling and logging patterns
- Always validate changes with appropriate tests before completion
