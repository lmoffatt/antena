cmake_minimum_required(VERSION 3.8)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
#set(CMAKE_CXX_STANDARD 17)

project(antena LANGUAGES CUDA CXX)


include_directories( quimulun/ )

set_source_files_properties(main.cpp PROPERTIES LANGUAGE CUDA)
#set_source_files_properties(main.cpp PROPERTIES LANGUAGE CXX)

add_executable(antena main.cpp  )

set(CMAKE_VERBOSE_MAKEFILE ON)

target_compile_options(antena PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-std=c++17 -g -ftemplate-backtrace-limit=0 -ftemplate-depth=1024 -expt-relaxed-constexpr >)

set_target_properties(antena PROPERTIES LINKER_LANGUAGE CXX)
