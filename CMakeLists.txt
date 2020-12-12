cmake_minimum_required(VERSION 3.17)
project(PIMC)

set(CMAKE_CXX_STANDARD 14)

add_executable(PIMC
        inc/harmonic_angle.h
        inc/harmonic_bond.h
        inc/coordinates.h
        inc/lennard_jones.h
        inc/physical_constants.h
        inc/electron.h
        inc/move_settings.h
        src/PIMC.cpp)

target_include_directories(PIMC PRIVATE inc)