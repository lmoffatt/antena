cmake_minimum_required(VERSION 3.5)
include_directories( quimulun/ )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
file(GLOB source_files "*.h" "*.cpp")

project(antena LANGUAGES CXX)
add_compile_options(-Wall -Wextra -Werror=return-type -std=c++17 -ftemplate-backtrace-limit=0 -fmax-errors=150 -lpthread -fopenmp)
add_executable(antena	 ${source_files})
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(antena PUBLIC OpenMP::OpenMP_CXX)
endif()
