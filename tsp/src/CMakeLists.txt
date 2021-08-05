cmake_minimum_required(VERSION 3.12)
project(tsp)
set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp Graph/Graph.cpp Graph/Graph.h UnionFind/UnionFind.h UnionFind/UnionFind.cpp Utils/utils.cpp)

add_executable(tsp ${SOURCE_FILES})