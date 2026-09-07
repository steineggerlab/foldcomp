# CudaBuild.cmake — all BUILD_CUDA-specific build configuration.
#
# Included from the top-level CMakeLists.txt only when BUILD_CUDA is ON, via
# include(CudaBuild) (resolved through CMAKE_MODULE_PATH).
#
# Ordering contract: this file MUST be included
#   * AFTER add_subdirectory(src)      — it appends to the foldcomp_* source lists, and
#   * BEFORE the foldcomp target is created — it enables the CUDA language and the target
#                                            compiles the appended .cu sources.
# include() runs in the includer's directory scope (no new variable scope), so the
# enable_language() call and the list(APPEND ...) below take effect at the top level.

# Record whether the user pinned CMAKE_CUDA_ARCHITECTURES BEFORE enable_language runs.
# Policy CMP0104 (NEW) auto-initializes CMAKE_CUDA_ARCHITECTURES to a default (e.g. "52")
# in the cache during enable_language — needed so configure-time try_compile()s (FindOpenMP
# etc.) work — which makes a later `if(NOT DEFINED ...)` check useless. So we decide here,
# while the cache is still clean, and remember it so reconfigures stay stable. A fresh
# configure with -DCMAKE_CUDA_ARCHITECTURES=... sets this ON; otherwise OFF and we discover.
if(NOT DEFINED FOLDCOMP_CUDA_ARCHS_USER_PINNED)
    if(DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(FOLDCOMP_CUDA_ARCHS_USER_PINNED ON CACHE INTERNAL "user pinned CMAKE_CUDA_ARCHITECTURES")
    else()
        set(FOLDCOMP_CUDA_ARCHS_USER_PINNED OFF CACHE INTERNAL "user pinned CMAKE_CUDA_ARCHITECTURES")
    endif()
endif()

# --- CUDA toolchain: must be at file/directory scope, never inside a function() ---
enable_language(CUDA)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# --- fold the CUDA sources into the core file lists (mutates includer-scope vars) ---
list(APPEND foldcomp_header_files ${foldcomp_cuda_header_files})
list(APPEND foldcomp_source_files ${foldcomp_cuda_source_files})

# ---------------------------------------------------------------------------
# CUDA architecture discovery
#
# Instead of hard-coding a fixed arch list (which breaks on toolkits that drop
# older archs — e.g. CUDA 13 removes sm_60/sm_70), probe the active compiler for
# the archs it actually accepts, then emit them in CMAKE_CUDA_ARCHITECTURES form.
# Mechanism ported from DALI's cmake/CUDA_utils.cmake
# (CUDA_check_cudacc_flag / CUDA_find_supported_arch_values / CUDA_get_cmake_cuda_archs).
# A user-supplied -DCMAKE_CUDA_ARCHITECTURES=... still wins and skips discovery.
# ---------------------------------------------------------------------------

# Check whether passing `flags` to the CUDA compiler succeeds. Unix only:
# equivalent to dry-running preprocessing on an empty .cu file and checking the
# exit code ($ nvcc ${flags} --dryrun -E -x cu /dev/null, or the clang variant).
function(CUDA_check_cudacc_flag out_status compiler flags)
    if(${compiler} MATCHES "clang")
        set(preprocess_empty_cu_file "-E" "-x" "cuda" "/dev/null")
    else()
        set(preprocess_empty_cu_file "--dryrun" "-E" "-x" "cu" "/dev/null")
    endif()
    execute_process(COMMAND ${compiler} ${flags} ${preprocess_empty_cu_file}
                    RESULT_VARIABLE tmp_out_status
                    OUTPUT_QUIET ERROR_QUIET)
    if(tmp_out_status EQUAL 0)
        set(${out_status} TRUE PARENT_SCOPE)
    else()
        set(${out_status} FALSE PARENT_SCOPE)
    endif()
endfunction()

# From a candidate arch list, keep only the values the given compiler accepts.
function(CUDA_find_supported_arch_values out_arch_values_allowed compiler arch_values_to_check)
    set(arch_list ${arch_values_to_check} ${ARGN})
    foreach(arch IN LISTS arch_list)
        if(${compiler} MATCHES "clang")
            CUDA_check_cudacc_flag(supported ${compiler} "--cuda-gpu-arch=sm_${arch}")
        else()
            CUDA_check_cudacc_flag(supported ${compiler} "-arch=sm_${arch}")
        endif()
        if(supported)
            set(out_list ${out_list} ${arch})
        endif()
    endforeach()
    set(${out_arch_values_allowed} ${out_list} PARENT_SCOPE)
endfunction()

# Turn a sorted arch list into CMAKE_CUDA_ARCHITECTURES form: "xx-real" for every
# arch, plus "<last>-virtual" so PTX is emitted for the newest arch, keeping
# forward compatibility with future GPUs.
function(CUDA_get_cmake_cuda_archs out_args_list arch_values)
    set(arch_list ${arch_values} ${ARGN})
    set(out "")
    foreach(arch IN LISTS arch_list)
        set(out "${out};${arch}-real")
    endforeach()
    list(GET arch_list -1 last_arch)
    set(out "${out};${last_arch}-virtual")
    set(${out_args_list} ${out} PARENT_SCOPE)
endfunction()

# Candidate archs to probe (overridable). Discovery narrows this to what the active
# toolkit supports, so it can list newer archs safely — unsupported ones are dropped.
set(FOLDCOMP_CUDA_TARGET_ARCHS "60;70;75;80;90;100;110;120"
    CACHE STRING "Candidate CUDA architectures probed against the active toolkit")

if(NOT FOLDCOMP_CUDA_ARCHS_USER_PINNED)
    CUDA_find_supported_arch_values(FOLDCOMP_CUDA_supported_archs
        ${CMAKE_CUDA_COMPILER} ${FOLDCOMP_CUDA_TARGET_ARCHS})
    if(NOT FOLDCOMP_CUDA_supported_archs)
        message(WARNING "No candidate CUDA arch was accepted by ${CMAKE_CUDA_COMPILER}; "
                        "falling back to sm_80. Set -DCMAKE_CUDA_ARCHITECTURES to override.")
        set(FOLDCOMP_CUDA_supported_archs "80")
    endif()
    message(STATUS "Foldcomp CUDA architectures (toolkit-supported): ${FOLDCOMP_CUDA_supported_archs}")
    CUDA_get_cmake_cuda_archs(FOLDCOMP_CUDA_ARCHITECTURES ${FOLDCOMP_CUDA_supported_archs})
endif()

# Configure the given target for CUDA: toolkit, defs, links, and architectures.
# The microtar/OpenMP/ZLIB dependencies are CUDA-gated only for the library/python
# builds; the executable build links them unconditionally (with its own flag-based
# OpenMP handling) so we must not touch them here.
function(foldcomp_configure_cuda_target target)
    find_package(CUDAToolkit REQUIRED)

    # Visibility for the CUDA include dirs and the dependency links: only the
    # installed library exposes an interface, so PUBLIC there, PRIVATE otherwise.
    if(BUILD_LIBRARY)
        set(_vis PUBLIC)
    else()
        set(_vis PRIVATE)
    endif()

    target_include_directories(${target} ${_vis} ${CUDAToolkit_INCLUDE_DIRS})
    target_compile_definitions(${target} PUBLIC FOLDCOMP_WITH_CUDA)
    target_link_libraries(${target} PRIVATE CUDA::cudart_static)
    if(NOT FOLDCOMP_CUDA_ARCHS_USER_PINNED)
        set_property(TARGET ${target} PROPERTY CUDA_ARCHITECTURES ${FOLDCOMP_CUDA_ARCHITECTURES})
    endif()

    if(BUILD_LIBRARY OR BUILD_PYTHON)
        include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib
                            ${CMAKE_CURRENT_SOURCE_DIR}/lib/microtar)
        if(NOT TARGET microtar)
            add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/microtar
                             ${CMAKE_CURRENT_BINARY_DIR}/lib/microtar)
        endif()
        target_link_libraries(${target} ${_vis} microtar)
        find_package(OpenMP REQUIRED)
        target_link_libraries(${target} ${_vis} OpenMP::OpenMP_CXX)
        target_compile_definitions(${target} PUBLIC OPENMP)
        find_package(ZLIB REQUIRED)
        target_link_libraries(${target} ${_vis} ZLIB::ZLIB)
        target_compile_definitions(${target} PRIVATE FOLDCOMP_WITH_ZLIB)
    endif()
endfunction()

# Standalone C++ decompression benchmark — counterpart to test/benchmark_cpu_gpu.py.
# Built only in the executable + CUDA config. Reuses the core (+ CUDA) sources like
# foldcomp_audit does, plus the structure reader so it can compress a PDB/CIF in-process.
function(foldcomp_add_cuda_bench)
    if(BUILD_LIBRARY OR BUILD_PYTHON OR EMSCRIPTEN)
        return()
    endif()
    find_package(CUDAToolkit REQUIRED)
    add_executable(foldcomp_bench
        ${foldcomp_header_files}
        ${foldcomp_source_files}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/structure_reader.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test/benchmark_cpu_gpu.cpp)
    target_include_directories(foldcomp_bench PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu
        ${CMAKE_CURRENT_SOURCE_DIR}/lib
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/gemmi
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/microtar
        ${CUDAToolkit_INCLUDE_DIRS})
    target_compile_definitions(foldcomp_bench PRIVATE
        FOLDCOMP_WITH_CUDA
        FOLDCOMP_WITH_MMCIF_OUTPUT
        FOLDCOMP_WITH_STRUCTURE_READER
        FOLDCOMP_WITH_ZLIB
        OPENMP
        _USE_MATH_DEFINES=1)
    target_link_libraries(foldcomp_bench PRIVATE
        microtar ZLIB::ZLIB CUDA::cudart_static)
    # Link OpenMP via raw flags rather than the OpenMP::OpenMP_CXX imported target:
    # that target's INTERFACE_COMPILE_OPTIONS embeds a
    # $<$<COMPILE_LANGUAGE:CXX>:SHELL:...> generator expression which some
    # Makefiles-generator/compiler combinations (e.g. gcc 13 on aarch64) fail to
    # evaluate, leaking the literal "$<1:...>" text into the compile command.
    # Mirrors the flag-string approach the foldcomp executable target already uses
    # in CMakeLists.txt for the same reason.
    if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        target_link_libraries(foldcomp_bench PRIVATE OpenMP::OpenMP_CXX)
    else()
        target_link_libraries(foldcomp_bench PRIVATE "${OpenMP_CXX_FLAGS}")
        target_compile_options(foldcomp_bench PRIVATE "${OpenMP_CXX_FLAGS}")
    endif()
    if(NOT FOLDCOMP_CUDA_ARCHS_USER_PINNED)
        set_property(TARGET foldcomp_bench PROPERTY CUDA_ARCHITECTURES ${FOLDCOMP_CUDA_ARCHITECTURES})
    endif()
endfunction()
