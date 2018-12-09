cmake_minimum_required(VERSION 3.1.0)
project(Chestnut CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-but-set-variable -Wno-reorder -Wno-maybe-uninitialized")

#Add all *.cpp and *.h files with GLOB_RECURSE:
file(GLOB_RECURSE SOURCES "src/*.cc")
file(GLOB_RECURSE HEADERS "src/*.h")

set(INCLUDE_DIRS "")
foreach (headerfile ${HEADERS})
	get_filename_component(dir ${headerfile} PATH)
	list (APPEND INCLUDE_DIRS ${dir})
endforeach()
list(REMOVE_DUPLICATES INCLUDE_DIRS)

#Default build in release mode
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "set build type" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG "-O2 -g")
set(CMAKE_C_FLAGS_DEBUG "-O2 -g")

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIRS})

add_custom_target(
	format
	COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && bash -c 'find src/ -name "*.cpp" | xargs -I % clang-format-4.0 -i %'
	COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && bash -c 'find src/ -name "*.h" | xargs -I % clang-format-4.0 -i %'
)
