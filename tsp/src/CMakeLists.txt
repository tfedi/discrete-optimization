cmake_minimum_required(VERSION 3.12)
project(tsp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=gnu++11 -ggdb3 -g")

add_executable(tsp main.cpp)