OpenMS Build System Helpers
======

This folder contains all scripts used by the OpenMS build system. Please add
only build system relevant scripts. All external scripts should go into
separate directories (e.g., `cmake/modules` for `Find*` scripts).

## Packaging Scripts

- `generate_applications_component_plist.cmake` - Generates ApplicationsComponent.plist for macOS PKG packages, listing all GUI application bundles with their relocatable paths.
