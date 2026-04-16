cmake_minimum_required(VERSION 3.15)

# Include Guard
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")

# Prevent in-source builds
if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL ERROR "Source folder cannot be safely used as build folder!")
endif()

##########################################################
# oncvpsp
##########################################################
project(metapsp
    VERSION 1.0.2
    DESCRIPTION "Accurate pseudopotentials for metagga calculations"
    LANGUAGES Fortran)

##########################################################
# Options
##########################################################
option(METAPSP_WITH_PSML "Enable PSML output support" ON)

###########################################################
# Include CMake modules
###########################################################
include(cmake/CPM.cmake)

##########################################################
# OS specific settings
##########################################################
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    message(STATUS "Configuring for macOS")
    execute_process(COMMAND sw_vers -productVersion
        OUTPUT_VARIABLE CMAKE_OSX_DEPLOYMENT_TARGET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "  macOS deployment target: ${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()

##########################################################
# Build types
##########################################################
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -finit-real=snan")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -finit-real=snan -DNDEBUG")

##########################################################
# Define the paths for static libraries and executables
##########################################################
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries."
)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries."
)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    CACHE
    PATH "Single output directory for building all executables."
)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
    CACHE
    PATH "Single output directory for building all Fortran modules."
)

###########################################################
# Common variables
###########################################################
set(METAPSP_DEPENDENCY_LIBRARIES "")

###########################################################
# LAPACK
###########################################################
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    # On macOS, FindLAPACK often finds the Accelerate framework,
    # which doesn't provide the correct Fortran API out of the box.
    # So, we do the following:
    # 1) If the user specifies LAPACK_ROOT, try to find LAPACK there
    # 2) If they don't set LAPACK_ROOT or searching there fails
    #    a) and BLA_VENDOR is set to "Apple"
    #       use lepus2589's accelerate-lapack recipe via CPM to build the
    #       proper LAPACK Fortran interface to the Accelerate framework
    #    b) and BLA_VENDOR is _not_ set to "Apple"
    #       compile reference LAPACK via CPM
    if (DEFINED LAPACK_ROOT AND (BLA_VENDOR STREQUAL "Apple"))
        message(FATAL_ERROR "Cannot set LAPACK_ROOT and BLA_VENDOR='Apple' at the same time!")
    endif()
    # First, if the user specifies a LAPACK_ROOT, try to find LAPACK there.
    if (DEFINED LAPACK_ROOT AND (NOT BLA_VENDOR STREQUAL "Apple"))
        find_package(LAPACK HINTS ${LAPACK_ROOT})
        list(APPEND METAPSP_DEPENDENCY_LIBRARIES ${LAPACK_LIBRARIES})
    endif()
    # If they don't, or searching in LAPACK_ROOT failed, use
    # lepus2589's accelerate-lapack recipe via CPM to build the
    # proper LAPACK Fortran interface to the Accelerate framework.
    if (NOT LAPACK_FOUND)
        if (BLA_VENDOR STREQUAL "Apple")
            CPMAddPackage("gh:lepus2589/accelerate-lapack@1.6.0")
            set(LAPACK_FOUND TRUE)
            list(APPEND METAPSP_DEPENDENCY_LIBRARIES LAPACK::LAPACK)
        else()
            CPMAddPackage("gh:Reference-LAPACK/lapack@3.12.1")
            set(LAPACK_FOUND TRUE)
            list(APPEND METAPSP_DEPENDENCY_LIBRARIES lapack)
        endif()
    endif()
else()
    if (DEFINED LAPACK_ROOT)
        find_package(LAPACK HINTS ${LAPACK_ROOT})
    else()
        find_package(LAPACK)
    endif()
    if (LAPACK_FOUND)
        message(STATUS "Found LAPACK version ${LAPACK_VERSION}")
        message(STATUS "  LAPACK libraries: ${LAPACK_LIBRARIES}")
        list(APPEND METAPSP_DEPENDENCY_LIBRARIES ${LAPACK_LIBRARIES})
    else()
        CPMAddPackage("gh:Reference-LAPACK/lapack@v3.12.1")
        list(APPEND METAPSP_DEPENDENCY_LIBRARIES LAPACK::LAPACK)
    endif()
endif()

###########################################################
# libxmlf90 (via CPM)
###########################################################
if (METAPSP_WITH_PSML)
    CPMAddPackage("https://gitlab.com/siesta-project/libraries/xmlf90.git#1.6.3")
    list(APPEND METAPSP_DEPENDENCY_LIBRARIES xmlf90-lib)
endif()

###########################################################
# Subdirectories
###########################################################
add_subdirectory(src)

###########################################################
# Installation
###########################################################
include(GNUInstallDirs)
if (METAPSP_WITH_PSML)
    install(TARGETS metapsp oncvpspm
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
else()
    install(TARGETS metapsp oncvpspmx
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
endif()
