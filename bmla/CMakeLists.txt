cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(bmla)

include_directories(/usr/local/include)

AUX_SOURCE_DIRECTORY(src SRC_FILES) 

file (GLOB HEADER_FILES "include/bmla_bits/*.hpp")
file (GLOB HEADER_FILES_FWDS "include/bmla_bits/fwds/*.hpp")

add_library(bmla SHARED 
	${SRC_FILES}
	${HEADER_FILES}
	${HEADER_FILES_FWDS}
)

target_include_directories(bmla PRIVATE include/bmla_bits include/bmla_bits/fwds)

find_library(ARMADILLO_LIB armadillo HINTS /usr/local/lib/)
target_link_libraries(bmla PUBLIC ${ARMADILLO_LIB})

install(TARGETS bmla DESTINATION lib)

install(FILES 
	include/bmla 
	DESTINATION include)

install(DIRECTORY 
	include/bmla_bits 
	DESTINATION include)