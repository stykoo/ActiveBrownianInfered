cmake_minimum_required(VERSION 3.0)

PROJECT(ActiveBrownian)

# Two executables, one with visualization, the other without
set(EXECUTABLE_NAME "ActiveBrownian")
set(EXECUTABLE_NAME_NOVISU "ActiveBrownian_novisu")

# Build in release mode by default
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Debug)
endif()
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "-W -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "-Ofast -g -pg") # Profiling

# This is not a strict requirement
set(CMAKE_CXX_STANDARD 17)

# Packages needed
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(HDF5 REQUIRED COMPONENTS C CXX)
# If any of these is not found, only the novisu executable will be compiled
find_package(Threads REQUIRED)
find_package(SFML COMPONENTS system window graphics)
include_directories(${HDF5_INCLUDE_DIRS})

find_package(MKL)
include_directories(${MKL_INCLUDE_DIRS})

if(MKL_FOUND)
	set(EXECUTABLE_NAME "ActiveBrownian_MKL")
	set(EXECUTABLE_NAME_NOVISU "ActiveBrownian_MKL_novisu")
endif()

# Executable with visualization
if(DEFINED CMAKE_THREAD_LIBS_INIT AND ${SFML_FOUND})
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

	if(MKL_FOUND)
		target_compile_definitions(${EXECUTABLE_NAME} PRIVATE USE_MKL)
    	target_link_libraries(${EXECUTABLE_NAME} -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group pthread dl)
	endif()

	# Linking
	target_link_libraries(
		${EXECUTABLE_NAME}
		${HDF5_LIBRARIES}
		${CMAKE_THREAD_LIBS_INIT}
		${SFML_LIBRARIES}
	)

	message("Threads and SMFL were found: the executable with vizualization will be built.")
else()
	message("Threads or SMFL was not found: only the executable without vizualization will be built.")
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
if(MKL_FOUND)
	target_compile_definitions(${EXECUTABLE_NAME_NOVISU} PRIVATE USE_MKL)
	target_link_libraries(${EXECUTABLE_NAME_NOVISU} -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group pthread dl)
endif()


target_link_libraries(
	${EXECUTABLE_NAME_NOVISU}
	${HDF5_LIBRARIES}
)
