cmake_minimum_required(VERSION 3.0.0)

project(kallisto)

include(GNUInstallDirs)

if(NOT MAX_KMER_SIZE)
    set(MAX_KMER_SIZE "32")
endif()

set(DO_ENABLE_AVX2 "")
set(DO_ENABLE_COMPILATION_ARCH "")

if(ENABLE_AVX2 MATCHES "OFF")
    add_compile_definitions("ENABLE_AVX2=OFF")
    set(DO_ENABLE_AVX2 "-DENABLE_AVX2=OFF")
endif(ENABLE_AVX2 MATCHES "OFF")

if(COMPILATION_ARCH MATCHES "OFF")
    add_compile_definitions("COMPILATION_ARCH=OFF")
    set(DO_ENABLE_COMPILATION_ARCH "-DCOMPILATION_ARCH=OFF")
endif(COMPILATION_ARCH MATCHES "OFF")

add_compile_definitions("MAX_KMER_SIZE=${MAX_KMER_SIZE}")


option(USE_HDF5 "Compile with HDF5 support" OFF) #OFF by default
option(USE_BAM "Compile with HTSLIB support" OFF) # OFF by default

if(USE_HDF5)
    cmake_policy(SET CMP0074 NEW)
    add_compile_definitions("USE_HDF5=ON")
endif(USE_HDF5)

if(NOT USE_BAM)
    add_compile_definitions("NO_HTSLIB=ON")
endif()

set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(CMAKE_CXX_FLAGS_PROFILE "-g")

# Set Release type for builds where CMAKE_BUILD_TYPE is unset
# This is usually a good default as this implictly enables
#
#   CXXFLAGS = -O3 -DNDEBUG
#
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(${CMAKE_VERSION} VERSION_LESS 3.1)
    # CMake versions before 3.1 do not understand CMAKE_CXX_STANDARD
    # remove this block once CMake >=3.1 has fixated in the ecosystem
    add_compile_options(-std=c++11)
else()
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-std=c++17 COMPILER_SUPPORTS_CXX17)
    if(COMPILER_SUPPORTS_CXX17)
        set(CMAKE_CXX_STANDARD 17)
    else()
        set(CMAKE_CXX_STANDARD 11)
    endif()

    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

#add_compile_options(-Wall -Wno-unused-function)
add_compile_options(-Wno-deprecated-declarations -Wno-subobject-linkage) # Suppress bifrost warning
set(PROJECT_BIFROST_CMAKE_CXX_FLAGS "-Wno-subobject-linkage -Wno-return-type") # Suppress bifrost warning

if(LINK MATCHES static)
    message("static build")
ELSE(LINK MATCHES shared)
    message("shared build")
ENDIF(LINK MATCHES static)


include(ExternalProject)
if (USE_BAM)
ExternalProject_Add(htslib
    PREFIX ${PROJECT_SOURCE_DIR}/ext/htslib
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/htslib
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND "" #autoreconf -i && autoheader &&  autoconf && ${PROJECT_SOURCE_DIR}/ext/htslib/configure
        #--prefix=${PREFIX} --disable-bz2 --disable-lzma --disable-libcurl
    BUILD_COMMAND "" #make lib-static
    INSTALL_COMMAND ""
)
endif(USE_BAM)

ExternalProject_Add(bifrost
    PREFIX ${PROJECT_SOURCE_DIR}/ext/bifrost
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/bifrost
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND mkdir -p build && cd build && cmake .. -DMAX_KMER_SIZE=${MAX_KMER_SIZE} -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_CXX_FLAGS=${PROJECT_BIFROST_CMAKE_CXX_FLAGS} ${DO_ENABLE_AVX2} ${DO_ENABLE_COMPILATION_ARCH}
    BUILD_COMMAND cd build && make
    INSTALL_COMMAND ""
)

if (ZLIBNG)
    message("zlib-ng enabled.")
    ExternalProject_Add(zlib-ng
    PREFIX ${PROJECT_SOURCE_DIR}/ext/zlib-ng
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/zlib-ng
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND mkdir -p zlib-ng && cd zlib-ng && cmake .. -DZLIB_COMPAT=ON -DZLIB_ENABLE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=${PREFIX}
    BUILD_COMMAND cd zlib-ng && make
    INSTALL_COMMAND ""
    )
endif(ZLIBNG)

if (USE_BAM)
    include_directories(${htslib_PREFIX}/src/htslib)
endif(USE_BAM)

include_directories(${EXT_PROJECTS_DIR}/bifrost/build/src)

ExternalProject_Get_Property(bifrost install_dir)
include_directories(${install_dir}/src)




add_subdirectory(src)
include_directories(${EXT_PROJECTS_DIR})

option(BUILD_TESTING "Build unit tests." OFF)
include(CTest)

if (BUILD_TESTING)
    add_subdirectory(${EXT_PROJECTS_DIR}/catch)

    # Includes Catch in the project:
    include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})

    add_subdirectory(unit_tests)
endif(BUILD_TESTING)

option(BUILD_FUNCTESTING "Build functional tests." OFF)

if (BUILD_FUNCTESTING)
    add_subdirectory(func_tests)
    message("Functional testing enabled.")
    add_custom_target(test
    COMMAND /bin/bash ./func_tests/runtests.sh
    DEPENDS ./src/kallisto
    )
endif(BUILD_FUNCTESTING)

# enable_testing()
# add_test(MainTest test/tests)
