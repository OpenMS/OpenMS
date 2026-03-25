# Release Process Analysis

## 1. Executive Summary
The OpenMS release process currently relies heavily on manual interventions, fragmented dependency management, and disconnected CI/CD workflows. This document analyzes the existing workflow to identify core pain points, aiming to establish a foundation for a fully automated, idempotent, and developer-friendly release pipeline.

## 2. Current Release Workflow
The typical release workflow involves several coordinated but disconnected steps:
1. **Preparation**: Updating version numbers across `CMakeLists.txt`, headers, and documentation.
2. **Dependency Alignment**: Ensuring `contrib` and `THIRDPARTY` packages are up-to-date and compatible across supported platforms.
3. **Build & Test Validation**: Relying on nightly Jenkins builds and GitHub Actions to verify stability on Linux, macOS, and Windows.
4. **Artifact Generation**: Building binaries, installers, and Python wheels manually or via semi-automated scripts.
5. **Distribution**: Uploading artifacts to package managers (Homebrew, Chocolatey, apt) and PyPI, often requiring manual approval or separate repository updates.
6. **Documentation Deployment**: Generating Doxygen and Sphinx documentation and manually deploying it to hosting platforms.

## 3. Repositories Involved
| Repository | Role |
| :--- | :--- |
| `OpenMS/OpenMS` | Main C++ library, TOPP tools, and GUI components. |
| `OpenMS/contrib` | External dependencies and build recipes. |
| `OpenMS/THIRDPARTY` | Vendored code and additional third-party binaries. |
| `OpenMS/OpenMS-release-scripts` | Custom scripts used for packaging and deployment. |
| Homebrew/Chocolatey forks | Package manager specific release formulas and manifests. |

## 4. Artifacts Generated
| Artifact Type | Platform | Distribution Channel |
| :--- | :--- | :--- |
| **Command Line Tools (TOPP)** | Linux, macOS, Windows | apt, Homebrew, Chocolatey, ZIP archives |
| **GUI Installers** | macOS, Windows | `.dmg` and `.exe` installers via GitHub Releases |
| **Python Bindings (`pyOpenMS`)** | Linux, macOS, Windows | PyPI (wheels for multiple Python versions) |
| **Docker Images** | Linux | GitHub Container Registry (GHCR) |
| **Documentation** | All | ReadTheDocs, OpenMS Website |

## 5. CI/CD Workflow Breakdown
* **GitHub Actions**: Primarily handles pull request validation, linting (`cpplint`), basic configuration checks, and some pyOpenMS building. Heavily impacted by runner updates and cache invalidation.
* **Jenkins**: Manages the nightlies, matrix testing across different compiler versions (GCC, Clang, MSVC), and heavy integration tests.
* **CDash**: Collects and visualizes test results from Jenkins and developer machines.
* **Release Branching**: Creating a `release/x.y.z` branch triggers specific workflows, but these are not always sufficient to build the final production-ready artifacts without intervention.

## 6. Pain Points
### Fragmentation
Dependency management is scattered across `apt`, Homebrew, Chocolatey, and the internal `contrib` / `THIRDPARTY` trees. This requires maintaining multiple build environments and replicating configuration logic.

### Manual Steps
Creating a release is not a single-button operation. Maintainers must manually bump versions, trigger specific Jenkins jobs, verify package distributions, and update external package manager repositories.

### Lack of Idempotency
Release scripts and UI builds often depend on external state (e.g., current network availability of dependencies). A failed release step often requires manual cleanup or a complete restart rather than safely resuming.

### Cache Inefficiency
CI caching is tightly coupled strictly to runner image versions. When GitHub updates a runner image, the build caches are invalidated, causing massive spikes in build times and resource consumption.

## 7. Developer Experience Issues
* **Steep Onboarding Curve**: New contributors find it difficult to test the full release pipeline locally.
* **Long Feedback Loops**: Because caching is inefficient, waiting for CI to validate packaging changes can take hours.
* **Opaque Failures**: Failures in the release pipeline often produce dense logs in Jenkins that are hard to decipher without deep institutional knowledge.

## 8. Key Observations
* The dual CI system (GitHub Actions + Jenkins) creates a split brain regarding build configuration and testing validation.
* `pyOpenMS` wrapping and wheel building have drastically different performance profiles and failure modes compared to the core C++ build.
* Fixing the release process will require decoupling the *build* of artifacts from their *distribution*.

## 9. Future Direction
The immediate next steps involve mapping out a unified CI/CD architecture based entirely on modern GitHub Actions, utilizing standardized metadata (like `pyproject.toml` and standard CMake presets). We should move toward an automated tag-triggered release pipeline that relies on strictly versioned compiler environments (via Docker) to ensure reproducible and idempotent artifact generation.
