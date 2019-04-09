cmake_minimum_required(VERSION 3.13)
project(MagneticSpheres)

set(CMAKE_CXX_STANDARD 14)

add_executable(MagneticSpheres main.cpp dipole.cpp dipole.h misc.cpp misc.h cluster.cpp cluster.h metropolis.cpp metropolis.h)