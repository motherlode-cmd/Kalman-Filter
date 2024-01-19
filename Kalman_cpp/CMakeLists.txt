cmake_minimum_required(VERSION 3.25)
project(Kalman_cpp)

set(CMAKE_CXX_STANDARD 17)

add_executable(Kalman_cpp main.cpp model/vector.hpp)
