# Header file-set regression tests

This standalone fixture uses the production library and installation helpers,
plus the actual test-framework target. It builds only a tiny probe library and
the standalone test framework, without OpenMS or its external dependencies.

Run with CMake 3.24 or newer and a C++23 compiler:

```sh
cmake -DBINARY_DIR=/tmp/openms-header-sets \
  -P src/tests/header_sets/run.cmake
```

Set `TEST_QT=ON` when Qt6 Core is available to check AUTOMOC on a header that
appears only in a public file set. `CONSUMER_CMAKE=/path/to/cmake-3.22` checks the
supported older consumer with real compilation and linking of both build-tree
and installed exports. `GENERATOR` defaults to Ninja; Unix Makefiles and
multi-configuration generators are supported. Compiler, toolchain and dependency
search settings can be passed as `-D` arguments. Logs are saved in `BINARY_DIR`.

The tests cover generated and source header contents, private-header exclusion,
default and explicit component installs, the test framework's `EXCLUDE_FROM_ALL`,
and the JSON guard with verification disabled. A deliberately broken public
header confirms that default builds skip verification and that explicitly
building `all_verify_interface_header_sets` detects an unavailable private
dependency. File sets do not enforce completeness of a component installation,
and compiler checks cannot reject private headers visible through shared or
system include directories; the configure-time JSON guard remains necessary.

`OPENMS_TEST_INSTALLED_CONSUMER` registers this fixture as the `HeaderSets` CTest
test. CI also runs the minimum-version producer/consumer combination separately.
The `linux-x64-ci` preset enables verification for the OpenMS-owned targets, and
the build action explicitly builds `all_verify_interface_header_sets`.
