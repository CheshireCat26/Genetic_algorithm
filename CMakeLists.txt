cmake_minimum_required(VERSION 3.15)
project(genetic_algorithm)

set(CMAKE_CXX_STANDARD 17)

add_executable(genetic_algorithm main.cpp Genetic_algorithm.cpp Genetic_algorithm.h)
add_subdirectory(Google_tests)