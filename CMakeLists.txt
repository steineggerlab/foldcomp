# Author: Milot Mirdita (milot@mirdita.de), Hyunbin Kim (khb7840@gmail.com)
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(foldcomp C CXX)

option(BUILD_LIBRARY "Build library" OFF)
option(BUILD_PYTHON "Build python support" OFF)
option(BUILD_WERROR "Build with warnings as errors" OFF)
option(BUILD_IWYU "Build with include-what-you-use" OFF)
option(GCS_SUPPORT "Enable Google Cloud Storage support" OFF)
option(ENABLE_AUDIT "Build foldcomp_audit helper" OFF)

include_directories(src)
add_subdirectory(src)

# For windows, include lib/windows
if(WIN32)
    include_directories(lib/windows)
endif(WIN32)

if(HAVE_SANITIZER)
    include(FindUBSan)
    include(FindASan)
    include(FindMSan)
    include(FindTSan)
endif()

if(BUILD_LIBRARY)
    add_library(foldcomp
        ${foldcomp_header_files}
        ${foldcomp_source_files})
elseif(BUILD_PYTHON)
        include_directories(lib)
        include_directories(lib/gemmi)
        if(MSVC)
            install(FILES lib/windows/dirent.h DESTINATION include)
        endif(MSVC)
        find_package(PythonInterp)
    	find_package(PythonLibs) # Trying to fix cibuildwheel ubuntu error
        find_package(PythonExtensions REQUIRED)
        add_library(foldcomp MODULE
            ${foldcomp_header_files}
            ${foldcomp_source_files}
            src/structure_reader.cpp
            foldcomp/foldcomp.cxx)
        target_compile_definitions(foldcomp PRIVATE FOLDCOMP_WITH_MMCIF_OUTPUT)
        target_compile_definitions(foldcomp PRIVATE FOLDCOMP_WITH_STRUCTURE_READER)
        python_extension_module(foldcomp)
        install(TARGETS foldcomp LIBRARY DESTINATION foldcomp)
else()
    include_directories(lib)

    add_executable(foldcomp
        ${foldcomp_header_files}
        ${foldcomp_source_files}
        ${foldcomp_exe_header_files}
        ${foldcomp_exe_source_files})

    target_compile_definitions(foldcomp PUBLIC FOLDCOMP_EXECUTABLE)
    target_compile_definitions(foldcomp PRIVATE FOLDCOMP_WITH_MMCIF_OUTPUT)
    target_compile_definitions(foldcomp PRIVATE FOLDCOMP_WITH_STRUCTURE_READER)
    # For debugging
    # target_compile_definitions(foldcomp PUBLIC _GLIBCXX_DEBUG=1 _LIBCPP_DEBUG=1)
    
    if(NOT EMSCRIPTEN)
        # For local compilation, openmp and zlib are required
        find_package(OpenMP REQUIRED)
        if(OPENMP_CXX_FOUND)
            if((CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
                target_link_libraries(foldcomp PUBLIC OpenMP::OpenMP_CXX)
            else()
                target_link_libraries(foldcomp PRIVATE "${OpenMP_CXX_FLAGS}")
                target_compile_options(foldcomp PRIVATE "${OpenMP_CXX_FLAGS}")
            endif()
            target_compile_definitions(foldcomp PUBLIC OPENMP)
        endif()
        find_package(ZLIB REQUIRED)
        target_link_libraries(foldcomp PUBLIC ZLIB::ZLIB)
        target_compile_definitions(foldcomp PRIVATE FOLDCOMP_WITH_ZLIB)
    else()
        # For webassembly, not using openmp. zlib is added to compile flags
        set_target_properties(
            foldcomp
            PROPERTIES
            COMPILE_FLAGS -sUSE_ZLIB=1
            LINK_FLAGS "-sUSE_ZLIB=1 -sEXPORTED_RUNTIME_METHODS=callMain,FS -sINVOKE_RUN=0 -sFILESYSTEM=1 -sALLOW_MEMORY_GROWTH=1 -sTOTAL_MEMORY=256MB -sENVIRONMENT=web -sMODULARIZE=1 -s EXPORT_ES6=1 -sEXPORT_NAME=createFoldcomp -sSINGLE_FILE=0 -sASSERTIONS=0")
    endif()

    include_directories(lib/gemmi)

    include_directories(lib/microtar)
    add_subdirectory(lib/microtar)
    target_link_libraries(foldcomp PUBLIC microtar)

    if(ENABLE_AUDIT)
        add_executable(foldcomp_audit
            ${foldcomp_header_files}
            ${foldcomp_source_files}
            src/structure_reader.cpp
            src/foldcomp_audit.cpp)
        target_link_libraries(foldcomp_audit PUBLIC ZLIB::ZLIB)
        target_compile_definitions(foldcomp_audit PRIVATE FOLDCOMP_WITH_MMCIF_OUTPUT)
        target_compile_definitions(foldcomp_audit PRIVATE FOLDCOMP_WITH_STRUCTURE_READER)
        target_compile_definitions(foldcomp_audit PRIVATE FOLDCOMP_WITH_ZLIB)
    endif()

    if(GCS_SUPPORT)
        find_package(google_cloud_cpp_storage REQUIRED)
        target_link_libraries(foldcomp PUBLIC google-cloud-cpp::storage)
        target_compile_definitions(foldcomp PUBLIC HAVE_GCS)
    endif()
endif()

# always enable math defines like M_PI
target_compile_definitions(foldcomp PUBLIC _USE_MATH_DEFINES=1)
if(TARGET foldcomp_audit)
    target_compile_definitions(foldcomp_audit PUBLIC _USE_MATH_DEFINES=1)
endif()

if(MSVC)
    target_compile_options(foldcomp PRIVATE /W4)
    if(TARGET foldcomp_audit)
        target_compile_options(foldcomp_audit PRIVATE /W4)
    endif()
else()
    target_compile_options(foldcomp PRIVATE -Wall -Wextra -Wpedantic)
    if(TARGET foldcomp_audit)
        target_compile_options(foldcomp_audit PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

if(BUILD_WERROR)
    if(MSVC)
        target_compile_options(foldcomp PRIVATE /WX)
        if(TARGET foldcomp_audit)
            target_compile_options(foldcomp_audit PRIVATE /WX)
        endif()
    else()
        target_compile_options(foldcomp PRIVATE -Werror)
        if(TARGET foldcomp_audit)
            target_compile_options(foldcomp_audit PRIVATE -Werror)
        endif()
    endif()
endif()

if(BUILD_IWYU)
    find_program(iwyu_path NAMES include-what-you-use iwyu REQUIRED)
    set_property(TARGET foldcomp PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path} -Xiwyu --mapping_file=${CMAKE_CURRENT_LIST_DIR}/cmake/include-what-you-use.imp)
endif()

enable_testing()
