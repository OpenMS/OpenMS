# GitHub Copilot Instructions for OpenMS

## Important Build Instructions

**DO NOT attempt to build OpenMS** - The build process is extremely resource-intensive and costly. Building OpenMS requires significant computational resources and time.

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

### Architecture Components

1. **OpenMS Library**
   - Core data structures (MSSpectrum, MSExperiment, Feature, FeatureMap)
   - Signal processing and analysis algorithms
   - File format handling and conversion

2. **TOPP Tools**
   - Command-line tools built on the OpenMS library
   - Standardized parameter handling via CTD scheme
   - Tools for identification, quantification, and data processing

3. **User Interfaces**
   - TOPPView for data visualization
   - TOPPAS for workflow creation
   - pyOpenMS for Python scripting

4. **pyOpenMS**
   - Python bindings generated via autowrap
   - Integration with scientific Python ecosystem
   - Located in `src/pyOpenMS/`

## Development Guidelines

1. **Code Style**
   - Follow the existing C++ coding conventions in the codebase
   - Use the established naming patterns for classes, methods, and variables
   - Maintain consistency with the surrounding code

2. **Testing**
   - Write unit tests for new functionality
   - Ensure existing tests pass before suggesting changes
   - Test files are located in `src/tests/`

3. **Documentation**
   - Add Doxygen comments for new public methods and classes
   - Update relevant documentation when modifying existing functionality

4. **Common Patterns**
   - Use OpenMS data structures (e.g., `MSExperiment`, `FeatureMap`, `PeptideIdentification`)
   - Follow the established error handling patterns
   - Utilize OpenMS logging mechanisms

5. **Performance Considerations**
   - Be mindful of memory usage when processing large datasets
   - Consider algorithmic complexity for data processing operations
   - Use appropriate OpenMS containers and algorithms

## Key Directories

- `src/openms/` - Core library source code
- `src/openms/include/OpenMS/` - Header files
- `src/topp/` - TOPP tools (command-line applications)
- `src/tests/` - Unit and integration tests
- `doc/` - Documentation

## Project Structure

```
OpenMS/
├── cmake/                   # CMake build system files
├── contrib/                 # Third-party dependencies
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
└── tools/                   # Development utilities
```

## Testing Infrastructure

- **Class tests**: Unit tests in `src/tests/class_tests/openms/`
- **TOPP tests**: Integration tests in `src/tests/topp/`
- **Python tests**: pyOpenMS tests in `src/pyOpenMS/tests/`
- Tests use the Google Test framework
- Follow naming convention: `ClassNameTest.cpp` for `ClassName`

## Assistance Focus

When providing code suggestions:
- Analyze existing code patterns before suggesting changes
- Respect the established architecture
- Consider cross-platform compatibility (Windows, Linux, macOS)
- Be aware of OpenMS dependencies and external libraries
- Use existing OpenMS data structures and algorithms
- Follow the established error handling patterns
- Maintain consistency with surrounding code style