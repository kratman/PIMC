cmake_minimum_required(VERSION 3.17)
project(PIMC)

set(CMAKE_CXX_STANDARD 14)

add_executable(PIMC
        inc/center_of_mass.h
        inc/coordinates.h
        inc/electron.h
        inc/harmonic_angle.h
        inc/harmonic_bond.h
        inc/lennard_jones.h
        inc/miscellaneous_constants.h
        inc/monte_carlo.h
        inc/move_settings.h
        inc/output.h
        inc/particle.h
        inc/path_integral.h
        inc/physical_constants.h
        inc/physics_settings.h
        inc/system.h
        inc/utilities.h
        src/PIMC.cpp)

target_include_directories(PIMC PRIVATE inc)