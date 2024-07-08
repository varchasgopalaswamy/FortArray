cmake_minimum_required(VERSION 3.28)
project(FortArray)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3  -fsave-optimization-record -march=native -Rpass-missed=.* -Rpass=.* -Rpass-analysis=.*")

add_executable(FortArray main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(FortArray PUBLIC OpenMP::OpenMP_CXX)
endif()
include_directories(include)