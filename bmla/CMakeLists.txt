cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(bmla)

add_library(bmla SHARED 
	src/connections.cpp
	src/general.cpp
	src/ixn_dicts.cpp
	src/ixn_param.cpp
	src/lattice.cpp
	src/moment.cpp
	src/opt_problem.cpp
	src/species.cpp
	src/unit.cpp
)

install(TARGETS bmla DESTINATION lib)

install(FILES 
	include/bmla 
	DESTINATION include)

install(DIRECTORY 
	include/bmla_bits 
	DESTINATION include)