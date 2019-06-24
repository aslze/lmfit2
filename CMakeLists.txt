cmake_minimum_required(VERSION 2.8.12)

project(lmfit2)

option(LMFIT_OPENMP "Use OpenMP" OFF)

if(LMFIT_OPENMP)
	find_package(OpenMP)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions(-DLMFIT_OPENMP)
endif()

set(TARGET lmfit2)

set(SRC
	main.cpp
	lmfit.h
	lmfit.c
)

add_executable(${TARGET} ${SRC})
