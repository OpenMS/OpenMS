project("OpenMS-seqan")
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# export include directory
set(SEQAN_INCLUDE_DIRECTORY
  ${PROJECT_SOURCE_DIR}
  CACHE
  INTERNAL "seqan include directories" FORCE
)
