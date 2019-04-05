cmake_minimum_required(VERSION 3.13)
project(MagneticSpheres)

set(CMAKE_CXX_STANDARD 14)

add_executable(MagneticSpheres main.cpp Dipole.cpp Dipole.h misc.cpp misc.h Cluster.cpp Cluster.h)