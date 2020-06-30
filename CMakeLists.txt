cmake_minimum_required(VERSION 3.8)
project(program)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp input.txt)
add_executable(program ${SOURCE_FILES})