# OpenMS Architecture Documentation

## Introduction

OpenMS is an open‐source software framework dedicated to mass spectrometry data analysis in proteomics, metabolomics, and related fields. It provides a robust, flexible, and extensible platform that supports rapid algorithm development, data processing, and the creation of custom analysis pipelines. This document outlines the core architectural components, design principles, and future directions of the OpenMS project.

## Architecture Principles

- **Modularity & Extensibility:** Built as a collection of loosely coupled libraries and tools, OpenMS allows developers to easily integrate new algorithms and extend existing functionality.
- **Performance & Scalability:** Written in modern C++ and optimized for parallel processing, OpenMS is designed to handle large datasets efficiently.
- **Portability:** OpenMS runs natively on Windows, Linux, and macOS.
- **Robustness & Reliability:** Extensive error handling, automated unit and integration tests, and continuous integration (e.g., via CDASH) ensure high-quality, reproducible analyses.
- **User-Centric Design:** With a suite of command-line tools (TOPP), graphical interfaces (TOPPView, TOPPAS), and Python bindings (pyOpenMS), OpenMS caters to both developers and end users.

## Components

### 1. OpenMS Library
- **Kernel Classes**: Essential data structures for representing mass spectrometry data
  - `Peak1D`: Basic data structure for mass-to-charge ratio (m/z) and intensity pairs
  - `ChromatogramPeak`: Data structure for retention time and intensity pairs
  - `MSSpectrum`: Represents a single mass spectrum
  - `MSExperiment`: Container for multiple spectra and chromatograms
  - `Feature`: Represents a detected feature in a mass spectrum, characterized by properties like m/z, retention time, and intensity.
  - `FeatureMap`: Container for features detected in a sample
  - `ConsensusMap`: Container for features identified across multiple samples or conditions
  - `PeptideIdentification`: Container for peptide identification results
  - `ProteinIdentification`: Container for protein identification results
  - ...

- **Algorithms:**  
  Implements key algorithms for signal processing, feature detection, quantification, peptide identification, protein inference, alignment, and others.
  
- **File Handling & Format Support:**  
  - Robust support for standard MS file formats (mzML, mzXML, mgf, mzIdentML, TraML, ...).
  - Utilities for file conversion and handling compressed data are also provided.

### 2. TOPP Tools (The OpenMS Pipeline Tools)

TOPP consists of comprehensive suite of command-line tools built on top of the OpenMS library that can be easily integrated into workflow systems.

- **Tool Architecture**:
  - Common interface for all tools, ensuring consistency.
  - Standardized parameter handling system for configuring tool behavior.
  - Logging and error reporting mechanisms for debugging and monitoring.
  - Progress monitoring to provide feedback during long-running processes.
  - Input/output file handling to manage data flow between tools.

### 3. User Interfaces
- **Graphical Tools:**  
  - **TOPPView:** A dedicated viewer for raw spectra, chromatograms, and identification results.
  - Raw data visualization
  - Spectrum and chromatogram viewing
  - 2D/3D data representation 
  - Identification result visualization
  - Interactive data analysis

### 4. Scripting  
- **Scripting & API: pyOpenMS** 
  - Python bindings that expose core functionalities
  - Rapid prototyping of algorithms in Python.
  -Integration with other scientific Python libraries (e.g., NumPy, SciPy).
  - Development of custom data processing workflows.
  - Interactive data analysis and visualization in Jupyter notebooks.

## Data flowchart

```mermaid	
flowchart TD
    A[Raw Data] --> B[Open file format]
    B --> C[Analysis Pipeline (e.g., set of TOPP Tools, pyOpenMS Scripts)]
    C --> D[Export]
    D --> E[Visualization and Downstream Processing]
```

## Extension Points

- **Algorithm Integration:**  
  Well-defined interfaces enable the addition of new data processing and analysis algorithms.
- **File Format Support:**  
  File handlers and extension hooks allow support for additional file formats and external code integration.
- **Tool Development:**  
  Developers can build new TOPP tools by subclassing common base classes and integrating with the standardized parameter and logging systems.
- **Workflow Customization:**  
  Users can combine OpenMS tools with custom scripts (e.g., via pyOpenMS).

## Build System

- **CMake Configuration:**  
  OpenMS uses a CMake-based build system that ensures platform-independent compilation and simplifies dependency management.
- **Automated Testing & CI:**  
  A comprehensive suite of unit tests, integration tests, and nightly builds (e.g., on CDASH) maintain code quality and facilitate rapid detection of issues.

## Parallel Processing

1. OpenMS uses OpenMP as parallelization backend

## Documentation Standards

- **Code Documentation:**  
  Doxygen-generated API documentation, inline comments, and consistent coding standards ensure clarity and maintainability.
- **User Documentation:**  
  User guides, tutorials, and example workflows.
- **Developer Guidelines:**  
  Clear contribution guidelines and coding conventions are provided to ensure consistency and quality across the codebase.

## Testing Strategy

- **Unit Testing:**  
  A comprehensive suite of unit tests covers the majority of core functionalities.
- **Integration Testing:**  
  Tests ensure that individual tools work together seamlessly and that workflows produce consistent results.

## Deployment and Distribution

OpenMS and its Python bindings (pyOpenMS) are distributed through several channels to suit different use cases and environments:

- **Standalone Installers:**  
  - For Windows and macOS, standalone installers (e.g., drag-and-drop installers for macOS) are provided for releases.

- **Bioconda:**  
  - OpenMS (and its library component `libopenms`) as well as its tools are available via the Bioconda channel.
  - The Python bindings, **pyOpenMS**, are available on Bioconda or pypi.

- **Container Images:**  
  - Docker and Singularity container images are provided through the OpenMS GitHub Container Registry as well as via BioContainers. These images bundle the OpenMS library, executables, and pyOpenMS so that users can deploy OpenMS in cloud or HPC environments with minimal setup.  


## Maintenance Guidelines

- **Code Reviews:**  
  All contributions undergo peer review to maintain quality and adherence to coding standards.
- **Release Management:**  
  Flexible release cycles with defined versioning protocols
- **Issue Tracking:**  
  Community-reported issues and feature requests are managed via GitHub issues


## Contributing

For detailed information about contributing to OpenMS, please refer to the CONTRIBUTING.md file in the repository.
