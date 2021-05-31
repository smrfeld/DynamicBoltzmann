cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(dblz)

include_directories(/usr/local/include)

AUX_SOURCE_DIRECTORY(src SRC_FILES) 

file (GLOB HEADER_FILES "include/dblz_bits/*.hpp")
file (GLOB HEADER_FILES_FWDS "include/dblz_bits/fwds/*.hpp")

add_library(bmla SHARED 
	src/center.cpp
	src/fname.cpp
	src/general.cpp
	src/ixn_param.cpp
	src/lattice_1d_fully_visible.cpp
	src/lattice_centered_hom.cpp
	src/lattice.cpp
	src/lattice_alternating_binary.cpp
	src/moment_diff.cpp
	src/opt_problem_static.cpp
	src/species.cpp
	include/dblz_bits/center.hpp
	include/dblz_bits/fname.hpp
	include/dblz_bits/general.hpp
	include/dblz_bits/ixn_param.hpp
	include/dblz_bits/lattice_1d_fully_visible.hpp
	include/dblz_bits/lattice_centered_hom.hpp
	include/dblz_bits/lattice.hpp
	include/dblz_bits/lattice_alternating_binary.hpp
	include/dblz_bits/moment_diff.hpp
	include/dblz_bits/opt_problem_static.hpp
	include/dblz_bits/species.hpp
	${HEADER_FILES_FWDS}
)

add_library(dblz SHARED 
	src/adjoint.cpp
	src/adjoint_params.cpp
	src/adjoint_params_centered_hom.cpp
	src/adjoint_obs.cpp
	src/center.cpp
	src/center_traj.cpp
	src/diff_eq_rhs.cpp
	src/fname_traj.cpp
	src/general.cpp
	src/ixn_param.cpp
	src/ixn_param_traj.cpp
	src/lattice_1d_fully_visible.cpp
	src/lattice_traj_1d_fully_visible.cpp
	src/lattice_centered_hom.cpp
	src/lattice.cpp
	src/lattice_alternating_binary.cpp
	src/lattice_traj.cpp
	src/lattice_traj_alternating_binary.cpp
	src/lattice_traj_centered_hom.cpp
	src/moment_diff.cpp
	src/opt_problem_dynamic.cpp
	src/species.cpp
	include/dblz_bits/adjoint.hpp
	include/dblz_bits/adjoint_params.hpp
	include/dblz_bits/adjoint_params_centered_hom.hpp
	include/dblz_bits/adjoint_obs.hpp
	include/dblz_bits/center.hpp
	include/dblz_bits/center_traj.hpp
	include/dblz_bits/diff_eq_rhs.hpp
	include/dblz_bits/fname_traj.hpp
	include/dblz_bits/general.hpp
	include/dblz_bits/ixn_param.hpp
	include/dblz_bits/ixn_param_traj.hpp
	include/dblz_bits/lattice_1d_fully_visible.hpp
	include/dblz_bits/lattice_traj_1d_fully_visible.hpp
	include/dblz_bits/lattice.hpp
	include/dblz_bits/lattice_alternating_binary.hpp
	include/dblz_bits/lattice_centered_hom.hpp
	include/dblz_bits/lattice_traj.hpp
	include/dblz_bits/lattice_traj_alternating_binary.hpp
	include/dblz_bits/lattice_traj_centered_hom.hpp
	include/dblz_bits/moment_diff.hpp
	include/dblz_bits/opt_problem_dynamic.hpp
	include/dblz_bits/species.hpp
	${HEADER_FILES_FWDS}
)

target_include_directories(bmla PRIVATE include/dblz_bits include/dblz_bits/fwds)
target_include_directories(dblz PRIVATE include/dblz_bits include/dblz_bits/fwds)

# Dependencies
include(cmake/CPM.cmake)
CPMAddPackage(
	URL "https://sourceforge.net/projects/arma/files/armadillo-9.200.8.tar.xz"
	NAME armadillo
	VERSION 9.200.8
	)
CPMAddPackage("gh:smrfeld/Q3-C1-Finite-Elements#3.0@3.0")

# Link
target_link_libraries(bmla PUBLIC armadillo)
target_link_libraries(dblz PUBLIC armadillo)
target_link_libraries(dblz PUBLIC q3c1)

install(TARGETS bmla DESTINATION lib)
install(TARGETS dblz DESTINATION lib)

install(FILES 
	include/bmla 
	DESTINATION include)

install(FILES 
	include/dblz 
	DESTINATION include)

install(DIRECTORY 
	include/dblz_bits 
	DESTINATION include)