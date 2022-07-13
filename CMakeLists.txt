# Author: Hyunbin Kim (khb7840@gmail.com)
# Reference: Martin's CMakeLists.txt & many websites

# 01. minimum version
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
# 02. project name & language
project(foldu_compression C CXX)

set(FRAMEWORK_ONLY 1 CACHE INTERNAL "" FORCE) # QUESTION: WHAT DOES THIS MEAN??

# 03. settings for google test
set(MY_GLOBAL_LIB_DIR /usr/local/lib)
set(MY_GLOBAL_INCLUDE /usr/local/include)
set(MY_INCLUDES ${PROJECT_SOURCE_DIR}/src)

include_directories(lib)
include_directories(${MY_GLOBAL_LIB_DIR})
include_directories(${MY_GLOBAL_INCLUDE})
include_directories(src)


# 03-1. find pthreads
find_package(Threads REQUIRED)
enable_testing()

# 04. add subdirectories
add_subdirectory(src)
# 05. add include directories

# 06. common variables
file(GLOB MY_SRCS ${PROJECT_SOURCE_DIR}/src/*.cpp)

# 07. enable AddressSanitizer
# 2022-05-27 22:43:37 Currently disabled due to error
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# 08. TARGET = MAIN
set(MAIN_SRC ${PROJECT_SOURCE_DIR}/main.cpp)
set(TARGET_NAME "foldcomp")
add_executable(${TARGET_NAME} ${MY_SRCS} ${MAIN_SRC})

# OPENMP
find_package(OpenMP QUIET)

if(OPENMP_FOUND)
    # For GCC we dont want to do this since it breaks macOS static builds
    # It will link libgomp.a internally (through -fopenmp I guess)
    # and also link libgomp.dylib thus breaking static builds
    if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_link_libraries(${TARGET_NAME} ${OpenMP_CXX_LIBRARIES})
    endif()

    append_target_property(${TARGET_NAME} COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    append_target_property(${TARGET_NAME} LINK_FLAGS ${OpenMP_CXX_FLAGS})
endif()

target_include_directories(${TARGET_NAME} PUBLIC
    ${MY_INCLUDES}
)