cmake_minimum_required(VERSION 3.16)
project(vep_annotator VERSION 1.0.0 LANGUAGES CXX)

# C++17 is required for std::optional and structured bindings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Find required packages
find_package(ZLIB REQUIRED)

# Find optional htslib for tabix support
option(WITH_HTSLIB "Enable tabix support via htslib for large VCF files" ON)

set(HAVE_HTSLIB FALSE)
if(WITH_HTSLIB)
    # Try pkg-config first
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(HTSLIB QUIET htslib)
    endif()

    # If pkg-config didn't find it, try manual search
    if(NOT HTSLIB_FOUND)
        find_path(HTSLIB_INCLUDE_DIR htslib/hts.h
            PATHS
                /usr/local/include
                /opt/homebrew/include
                /usr/include
        )
        find_library(HTSLIB_LIBRARY hts
            PATHS
                /usr/local/lib
                /opt/homebrew/lib
                /usr/lib
        )
        if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
            set(HTSLIB_FOUND TRUE)
            set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR})
            set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})
        endif()
    endif()

    if(HTSLIB_FOUND)
        set(HAVE_HTSLIB TRUE)
        message(STATUS "Found htslib: ${HTSLIB_LIBRARIES}")
    else()
        message(STATUS "htslib not found - tabix support disabled")
        message(STATUS "Install htslib for large VCF support: brew install htslib (macOS)")
    endif()
endif()

# Find optional libBigWig for conservation scores
option(WITH_BIGWIG "Enable bigWig support for conservation scores" ON)

set(HAVE_BIGWIG FALSE)
if(WITH_BIGWIG)
    find_path(BIGWIG_INCLUDE_DIR bigWig.h
        PATHS
            /usr/local/include
            /opt/homebrew/include
            /usr/include
    )
    find_library(BIGWIG_LIBRARY BigWig
        PATHS
            /usr/local/lib
            /opt/homebrew/lib
            /usr/lib
    )
    if(BIGWIG_INCLUDE_DIR AND BIGWIG_LIBRARY)
        set(HAVE_BIGWIG TRUE)
        set(BIGWIG_INCLUDE_DIRS ${BIGWIG_INCLUDE_DIR})
        set(BIGWIG_LIBRARIES ${BIGWIG_LIBRARY})
        message(STATUS "Found libBigWig: ${BIGWIG_LIBRARIES}")
    else()
        message(STATUS "libBigWig not found - conservation score support disabled")
        message(STATUS "Install libBigWig for PhyloP/PhastCons/GERP support")
    endif()
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${ZLIB_INCLUDE_DIRS})

# Library sources
set(LIB_SOURCES
    src/vep_annotator.cpp
    src/file_parsers.cpp
    src/plugin_loader.cpp
    src/sources/dbnsfp_source.cpp
    src/sources/spliceai_source.cpp
    src/sources/maxentscan_source.cpp
    src/sources/dbscsnv_source.cpp
    src/sources/conservation_source.cpp
    src/sources/regulatory_source.cpp
    src/sources/domain_source.cpp
    src/sources/lof_source.cpp
)

# Create library
add_library(vep_annotator_lib STATIC ${LIB_SOURCES})
target_link_libraries(vep_annotator_lib
    PUBLIC
        ${ZLIB_LIBRARIES}
        ${CMAKE_DL_LIBS}
)
target_include_directories(vep_annotator_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Add htslib if found
if(HAVE_HTSLIB)
    target_compile_definitions(vep_annotator_lib PUBLIC HAVE_HTSLIB)
    target_include_directories(vep_annotator_lib PUBLIC ${HTSLIB_INCLUDE_DIRS})
    target_link_libraries(vep_annotator_lib PUBLIC ${HTSLIB_LIBRARIES})
endif()

# Add libBigWig if found
if(HAVE_BIGWIG)
    target_compile_definitions(vep_annotator_lib PUBLIC HAVE_BIGWIG)
    target_include_directories(vep_annotator_lib PUBLIC ${BIGWIG_INCLUDE_DIRS})
    target_link_libraries(vep_annotator_lib PUBLIC ${BIGWIG_LIBRARIES})
endif()

# Create main executable
add_executable(vep_annotator src/main.cpp)
target_link_libraries(vep_annotator
    PRIVATE
        vep_annotator_lib
)

# Create filter_vep executable
add_executable(filter_vep src/filter_vep.cpp)
target_include_directories(filter_vep PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(filter_vep
    PRIVATE
        ${ZLIB_LIBRARIES}
)

# Installation
install(TARGETS vep_annotator filter_vep DESTINATION bin)
install(TARGETS vep_annotator_lib DESTINATION lib)
install(FILES
    include/vep_annotator.hpp
    include/annotation_source.hpp
    include/annotation_sources.hpp
    include/file_parsers.hpp
    include/plugin.hpp
    include/output_writer.hpp
    include/transcript_filter.hpp
    include/structural_variant.hpp
    include/hgvs_parser.hpp
    include/gene_constraint.hpp
    include/exon_intron_numbers.hpp
    include/filter_vep.hpp
    DESTINATION include)

# Print configuration summary
message(STATUS "")
message(STATUS "VEP Annotator Configuration Summary")
message(STATUS "====================================")
message(STATUS "C++ Standard:      C++${CMAKE_CXX_STANDARD}")
message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "ZLIB include:      ${ZLIB_INCLUDE_DIRS}")
message(STATUS "ZLIB libraries:    ${ZLIB_LIBRARIES}")
if(HAVE_HTSLIB)
    message(STATUS "htslib:            ENABLED (tabix support for large VCF files)")
    message(STATUS "htslib include:    ${HTSLIB_INCLUDE_DIRS}")
    message(STATUS "htslib libraries:  ${HTSLIB_LIBRARIES}")
else()
    message(STATUS "htslib:            DISABLED (--annotation-tabix will not work)")
endif()
if(HAVE_BIGWIG)
    message(STATUS "libBigWig:         ENABLED (conservation scores)")
    message(STATUS "libBigWig include: ${BIGWIG_INCLUDE_DIRS}")
    message(STATUS "libBigWig library: ${BIGWIG_LIBRARIES}")
else()
    message(STATUS "libBigWig:         DISABLED (PhyloP/PhastCons/GERP not available)")
endif()
message(STATUS "")
message(STATUS "This is a pure C++ local implementation.")
message(STATUS "No external API calls - requires local GTF and FASTA files.")
message(STATUS "")
