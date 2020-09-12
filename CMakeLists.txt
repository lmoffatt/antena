cmake_minimum_required(VERSION 3.8)
#set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

project(antena LANGUAGES CUDA CXX)


include_directories( quimulun/ /usr/local/cuda-11.0/targets/x86_64-linux/include/)

#set_source_files_properties(main.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(main.cpp PROPERTIES LANGUAGE CXX)

add_executable(antena main.cpp antena_data.h  )

set(CMAKE_VERBOSE_MAKEFILE ON)

add_compile_options(-Wall -Wextra -Werror=return-type -std=c++17 -ftemplate-backtrace-limit=0 -fmax-errors=20 -fopenmp -fconstexpr-depth=30 )

#target_compile_options(antena PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-std=c++17 -g -ftemplate-backtrace-limit=0 -ftemplate-depth=1024 -expt-relaxed-constexpr >)

set_target_properties(antena PROPERTIES LINKER_LANGUAGE CXX)
