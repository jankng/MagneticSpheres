cmake_minimum_required(VERSION 3.13)
project(MagneticSpheres)

set(CMAKE_CXX_STANDARD 14)

add_executable(MagneticSpheres main.cpp dipole.cpp dipole.h misc.cpp misc.h cluster.cpp cluster.h metropolis.cpp metropolis.h gsledits.cpp gsledits.h)


# link with gsl library
#set(GSL_PATH ../dependencies/gsl)
include_directories(${GSL_PATH}/include)
#link_directories(${GSL_PATH}/bin)
target_link_libraries(MagneticSpheres gsl gslcblas pthread)