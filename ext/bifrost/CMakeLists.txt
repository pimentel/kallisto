cmake_minimum_required(VERSION 3.0.0)

project(Bifrost)

find_package(Threads REQUIRED)

# To enable a larger default k-mer size, replace MAX_KMER_SIZE with a larger multiple of 32: actual maximum k-mer size will be MAX_KMER_SIZE-1.
SET(MAX_GMER_SIZE "${MAX_KMER_SIZE}" CACHE STRING "MAX_GMER_SIZE")
# Enable architecture optimizations
SET(COMPILATION_ARCH "native" CACHE STRING "COMPILATION_ARCH")
# Enable AVX2 instructions
SET(ENABLE_AVX2 "ON" CACHE STRING "ENABLE_AVX2")

# Set some default compile flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set_property(SOURCE BlockedBloomFilter.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -funroll-loops")


#check if we are on arm64 and apple, if so, disable AVX2
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm")
	message("Disabling AVX2 instructions on arm64")
	set(ENABLE_AVX2 "OFF")
	set(COMPILATION_ARCH "OFF")
endif(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm")

if(COMPILATION_ARCH MATCHES "OFF")
	message("Disabling native architecture compilation (including AVX2)")
else(COMPILATION_ARCH MATCHES "OFF")
	message("Compilation architecture: ${COMPILATION_ARCH}")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=${COMPILATION_ARCH}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${COMPILATION_ARCH}")
endif(COMPILATION_ARCH MATCHES "OFF")

if(ENABLE_AVX2 MATCHES "OFF")
	message("Disabling AVX2 instructions")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mno-avx2")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-avx2")
endif(ENABLE_AVX2 MATCHES "OFF")

# Manages build types
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
	message("Build type: Debug")
	add_compile_options(-g)
else(CMAKE_BUILD_TYPE MATCHES Debug)
	if(CMAKE_BUILD_TYPE MATCHES Profile)
		message("Build type: Profiling")
		add_compile_options(-pg)
		set(CMAKE_SHARED_LINKER_FLAGS "-pg")
		set(CMAKE_EXE_LINKER_FLAGS "-pg")
	else(CMAKE_BUILD_TYPE MATCHES Profile)
		message("Build type: Release")
		add_compile_options(-O3)
	endif(CMAKE_BUILD_TYPE MATCHES Profile) 
endif(CMAKE_BUILD_TYPE MATCHES Debug)

MATH(EXPR PRINT_MAX_KMER_SIZE "${MAX_KMER_SIZE}-1")
message("Maximum k-mer size: " ${PRINT_MAX_KMER_SIZE})

MATH(EXPR PRINT_MAX_GMER_SIZE "${MAX_GMER_SIZE}-1")
message("Maximum g-mer size: " ${PRINT_MAX_GMER_SIZE})

add_subdirectory(src)
