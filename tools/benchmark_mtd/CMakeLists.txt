cmake_minimum_required(VERSION 3.24)
project(benchmark_mtd LANGUAGES CXX)

find_package(OpenMS REQUIRED CONFIG)
find_package(OpenMP REQUIRED COMPONENTS CXX)

add_executable(benchmark_mtd benchmark_mtd.cpp)
target_compile_features(benchmark_mtd PRIVATE cxx_std_23)
target_link_libraries(benchmark_mtd PRIVATE OpenMS::OpenMS OpenMP::OpenMP_CXX)
