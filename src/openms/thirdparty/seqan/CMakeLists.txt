project("OpenMS-seqan")
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

# export include directory
set(SEQAN_INCLUDE_DIRECTORY
  ${PROJECT_SOURCE_DIR}/include
  CACHE
  INTERNAL "seqan include directories" FORCE
)
