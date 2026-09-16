# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause

### collect the sources and headers of the OpenMS_CLI library
set(OpenMS_CLI_sources CACHE INTERNAL "This variable should hold all OpenMS_CLI sources at the end of the config step")
include(source/APPLICATIONS/sources.cmake)

set(OpenMS_CLI_sources_h CACHE INTERNAL "This variable should hold all OpenMS_CLI headers at the end of the config step")
include(include/OpenMS/APPLICATIONS/sources.cmake)
