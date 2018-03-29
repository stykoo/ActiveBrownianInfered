cmake_minimum_required(VERSION 2.8)

PROJECT(ActiveBrownian)

# Two executables, one with visualization, the other without
set(EXECUTABLE_NAME "ActiveBrownian")
set(EXECUTABLE_NAME_NOVISU "ActiveBrownian_novisu")

# Build in release mode by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
if (NOT EXECUTABLE_OUTPUT_PATH)
	if(CMAKE_BUILD_TYPE STREQUAL Release)
		set(EXECUTABLE_OUTPUT_PATH bin/)
	else()
		set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
	endif()
endif()

set(CMAKE_CXX_FLAGS "-W -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "-Ofast -g -pg") # Profiling

# This is not a strict requirement
set(CMAKE_CXX_STANDARD 14)


# Packages needed
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(Boost REQUIRED COMPONENTS program_options)
# If any of these is not found, only the novisu executable will be compiled
find_package(Threads REQUIRED)
find_package(SFML COMPONENTS system window graphics)
find_package(VTK)
include_directories(${HDF5_INCLUDE_DIRS})
include(${VTK_USE_FILE})

# Executable with visualization
if(DEFINED CMAKE_THREAD_LIBS_INIT AND ${SFML_FOUND} AND ${VTK_FOUND})
	file(
		GLOB
		source_files
		src/*
		src/visu/*
	)

	add_executable(
		${EXECUTABLE_NAME}
		${source_files}
	)

	# Linking
	target_link_libraries(
		${EXECUTABLE_NAME}
		${HDF5_LIBRARIES}
		${Boost_LIBRARIES}
		${CMAKE_THREAD_LIBS_INIT}
		${SFML_LIBRARIES}
		${VTK_LIBRARIES}
	)

	message("Threads, SMFL and VTK were found: the executable with vizualization will be built.")
else()
	message("Threads, SMFL or VTK was not found: only the executable without vizualization will be built.")
endif()


# Executable without visualization
file(
	GLOB
	source_files_novisu
	src/*
)

add_executable(
	${EXECUTABLE_NAME_NOVISU}
	${source_files_novisu}
)

# defines NOVISU
target_compile_definitions(${EXECUTABLE_NAME_NOVISU} PRIVATE NOVISU)

target_link_libraries(
	${EXECUTABLE_NAME_NOVISU}
	${HDF5_LIBRARIES}
	${Boost_LIBRARIES}
)
