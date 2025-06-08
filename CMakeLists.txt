cmake_minimum_required(VERSION 3.31)
project(untitled2)

set(CMAKE_CXX_STANDARD 20)

add_executable(untitled2 main.cpp
        params.h
        utils.h
        utils.cpp
        analytical.cpp
        analytical.h
        CALERF.cpp
        CALERF.h
        explicit.cpp
        explicit.h
        implicit.cpp
        implicit.h
        linalg.cpp
        linalg.h
        output.cpp
        output.h
        solver.cpp
        solver.h)
