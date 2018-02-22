cmake_minimum_required(VERSION 2.8)

PROJECT(ActiveBrownian)
set(EXECUTABLE_NAME "ActiveBrownian")
set(EXECUTABLE_OUTPUT_PATH bin/)

# Build in release mode by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-W -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# This is not a strict requirement
set(CMAKE_CXX_STANDARD 14)

# Packages needed
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(Threads)
find_package(Boost REQUIRED program_options)
find_package(SFML REQUIRED system window graphics)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Configuration of the executable
file(
	GLOB_RECURSE
	source_files
	src/*
)
add_executable(
	${EXECUTABLE_NAME}
	${source_files}
)

# Linking
target_link_libraries(
	${EXECUTABLE_NAME}
	${CMAKE_THREAD_LIBS_INIT}
	${Boost_LIBRARIES}
	${SFML_LIBRARIES}
	${VTK_LIBRARIES}
)
