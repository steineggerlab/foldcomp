# Author: Hyunbin Kim (khb7840@gmail.com)
# Reference: MMSeq2's CMakeLists.txt & many websites

# 01. minimum version
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# 02. project name & language
project(foldu_compression C CXX)

set(HAVE_GCS 0 CACHE BOOL "Have Google Cloud Storage SDK")

# 03. settings for google test
set(MY_GLOBAL_LIB_DIR /usr/local/lib)
set(MY_GLOBAL_INCLUDE /usr/local/include)
set(MY_INCLUDES ${PROJECT_SOURCE_DIR}/src)

include_directories(lib)
include_directories(${MY_GLOBAL_LIB_DIR})
include_directories(${MY_GLOBAL_INCLUDE})
include_directories(src)

# 03-1. find pthreads
# find_package(Threads REQUIRED)
enable_testing()

# 04. add subdirectories
add_subdirectory(src)
# 05. add include directories

# 06. common variables
file(GLOB MY_SRCS ${PROJECT_SOURCE_DIR}/src/*.cpp)

# 07. enable AddressSanitizer
if (HAVE_SANITIZER)
    include(FindUBSan)
    include(FindASan)
    include(FindMSan)
    include(FindTSan)
endif ()

# 08. TARGET = MAIN
set(MAIN_SRC ${PROJECT_SOURCE_DIR}/main.cpp)
set(TARGET_NAME "foldcomp")
add_executable(${TARGET_NAME} ${MY_SRCS} ${MAIN_SRC})

# 09. OPENMP
find_package(OpenMP REQUIRED)
if(OPENMP_CXX_FOUND)
    target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(${TARGET_NAME} PUBLIC
    ${MY_INCLUDES}
)

# 10. Gemmi - for reading cif files & gzipped files
find_package(ZLIB REQUIRED)
target_link_libraries(${TARGET_NAME} PUBLIC ZLIB::ZLIB)
include_directories(lib/gemmi)

# 11. Microtar - for reading/writing tar files
include_directories(lib/microtar)
add_subdirectory(lib/microtar)
target_link_libraries(${TARGET_NAME} PUBLIC microtar)

# 12
if(HAVE_GCS)
    find_package(google_cloud_cpp_storage REQUIRED)
    target_link_libraries(${TARGET_NAME} PUBLIC google-cloud-cpp::storage)
    target_compile_definitions(${TARGET_NAME} PUBLIC HAVE_GCS)
endif()
