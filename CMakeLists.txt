cmake_minimum_required(VERSION 3.25)
project(project)

set(CMAKE_CXX_STANDARD 17)

add_executable(project main.cpp)
add_definitions(-D_UNICODE -DUNICODE)