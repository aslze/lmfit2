cmake_minimum_required(VERSION 2.8.12)

project(lmfit2)

option(LMFIT_OPENMP "Use OpenMP")

if(LMFIT_OPENMP)
	find_package(OpenMP)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions(-DLMFIT_OPENMP)
endif()

set(TARGET lmfit2)

add_library(${TARGET} STATIC ${SRC})
target_include_directories(${TARGET} PUBLIC .)

set(SRC
	lmfit.h
	lmfit.c
)


if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
	
	set(TARGET lmfit-demo)

	set(SRC
		main.cpp
	)

	add_executable(${TARGET} ${SRC})
	target_link_libraries(${TARGET} lmfit2)

endif()

