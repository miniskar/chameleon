###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  @project CHAMELEON
#  CHAMELEON is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
# @version 0.9.2
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2014-11-16
#
###
cmake_minimum_required(VERSION 3.3)
project(CHAMELEON Fortran C CXX)

# directly make an error if in-source build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.\n"
   "Please create a build directory first and execute cmake configuration from "
   "this directory. Example: mkdir build && cd build && cmake ..")
endif()

# set project version number
set(CHAMELEON_VERSION_MAJOR 0)
set(CHAMELEON_VERSION_MINOR 9)
set(CHAMELEON_VERSION_MICRO 2)

set(CHAMELEON_CMAKE_DIR "" CACHE PATH "Directory of CHAMELEON CMake modules, can be external to the project")

# Add extra cmake module path and initialize morse cmake modules
# --------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
if(CHAMELEON_CMAKE_DIR)
    list(APPEND CMAKE_MODULE_PATH "${CHAMELEON_CMAKE_DIR}/cmake_modules/morse_cmake/modules")
    set(CHAMELEON_CMAKE_MODULE_PATH ${CHAMELEON_CMAKE_DIR}/cmake_modules/morse_cmake/modules )
elseif(EXISTS "${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake")
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules")
    set(CHAMELEON_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules )
else()
  message(FATAL_ERROR "Chameleon CMake system relies on morse_cmake modules developed here: "
    "https://gitlab.inria.fr/solverstack/morse_cmake. Please set CHAMELEON_CMAKE_DIR to this source "
    "directory.")
endif()

include(MorseInit)
include(GenPkgConfig)

# Ensure that the fortran compiler and c compiler specified are compatible
# ------------------------------------------------------------------------
include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/include/chameleon/mangling.h
                         MACRO_NAMESPACE "CHAMELEON_")

#############################################
#                                           #
#        Compilation of CHAMELEON           #
#                                           #
#############################################


###############################################################################
# Parameters/Options #
######################

# Set the RPATH config
# --------------------

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Misc options
# ------------
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Define precision supported by CHAMELEON
# -----------------------------------------
set( RP_CHAMELEON_DICTIONNARY ${CMAKE_SOURCE_DIR}/cmake_modules/local_subs.py )
set( RP_CHAMELEON_PRECISIONS  "s;d;c;z" )
include(RulesPrecisions)

option(CHAMELEON_USE_MIGRATE
  "This options enables the data migration in QR algorithms" ON)

# Options to select the runtime
# -----------------------------

# Create a list of possible runtime
set(CHAMELEON_SCHED_list PARSEC STARPU QUARK OPENMP
  CACHE INTERNAL "List of available runtimes" )

set( CHAMELEON_SCHED_PARSEC OFF CACHE INTERNAL
  "Enable PaRSEC scheduler as the default runtime
    (Conflict with other CHAMELEON_SCHED_* options)")
set(CHAMELEON_SCHED_STARPU OFF CACHE INTERNAL
    "Enable StarPU scheduler as the default runtime
    (Conflict with other CHAMELEON_SCHED_* options)")
set(CHAMELEON_SCHED_QUARK OFF CACHE INTERNAL
  "Enable Quark scheduler as the default runtime
    (Conflict with other CHAMELEON_SCHED_* options)")
set(CHAMELEON_SCHED_OPENMP OFF CACHE INTERNAL
  "Enable OpenMP scheduler as the default runtime
    (Conflict with other CHAMELEON_SCHED_* options)")

set( CHAMELEON_SCHED STARPU
  CACHE STRING "Choose the chameleon internal runtime from ${CHAMELEON_SCHED_list}")
set_property(CACHE CHAMELEON_SCHED PROPERTY STRINGS ${CHAMELEON_SCHED_list})

if ( NOT ${CHAMELEON_SCHED} IN_LIST CHAMELEON_SCHED_list )
  message(FATAL_ERROR "CHAMELEON_SCHED must be one of ${CHAMELEON_SCHED_list} (${CHAMELEON_SCHED} is incorrect)" )
endif()

set( CHAMELEON_SCHED_${CHAMELEON_SCHED} ON )

# Set default to StarPU if nothing specific is required by the user
message("-- ${BoldGreen}CHAMELEON_SCHED is set to use ${CHAMELEON_SCHED} runtime" )

# Check that one, and only one, SCHED option is set to ON
# count number of runtime sets to ON
math(EXPR number_of_active_runtime 0)
foreach (runtime ${CHAMELEON_SCHED_list})
    if (CHAMELEON_SCHED_${runtime})
        math(EXPR number_of_active_runtime "${number_of_active_runtime}+1")
    endif()
endforeach()
if (NOT number_of_active_runtime STREQUAL 1)
    message(FATAL_ERROR
            "Number of active runtime is ${number_of_active_runtime}, "
            "the user should activate one (and only one) runtime. ")
endif()

# Use intermediate variable since cmake_dependent_option doesn't have OR conditions
set(CHAMELEON_ENABLE_MPI OFF CACHE INTERNAL "Tells if MPI might be supported by the runtime")
if ( CHAMELEON_SCHED_PARSEC OR CHAMELEON_SCHED_STARPU )
     set(CHAMELEON_ENABLE_MPI ON FORCE)
endif()

# Use intermediate variable since cmake_dependent_option doesn't have OR conditions
set(CHAMELEON_ENABLE_CUDA OFF CACHE INTERNAL "Tells if CUDA might be supported by the runtime")
if ( CHAMELEON_SCHED_PARSEC OR CHAMELEON_SCHED_STARPU )
     set(CHAMELEON_ENABLE_CUDA ON FORCE)
endif()

# Additional options
# ------------------

# Enable the distributed interface (allowed only when StarPU or PaRSEC is enabled)
# TODO: Default should be changed to ON/OFF when it will be ok
cmake_dependent_option(CHAMELEON_USE_MPI
                    "Enable distributed memory through MPI" OFF
                    "CHAMELEON_ENABLE_MPI" OFF)
if (CHAMELEON_ENABLE_MPI AND NOT CHAMELEON_USE_MPI)
    message("-- ${BoldGreen}CHAMELEON_USE_MPI is set to OFF, turn it ON to use MPI (unsupported by Quark)${ColourReset}")
endif()

# Enable Cuda kernels if StarPU (only if StarPU or PaRSEC runtime is enabled)
# TODO: Default should be changed to ON/OFF when it will be ok
cmake_dependent_option(CHAMELEON_USE_CUDA
                    "Enable CUDA kernels" OFF
                    "CHAMELEON_ENABLE_CUDA" OFF)
if (CHAMELEON_ENABLE_CUDA AND NOT CHAMELEON_USE_CUDA)
    message("-- ${BoldGreen}CHAMELEON_USE_CUDA is set to OFF, turn it ON to use CUDA (unsupported by Quark)${ColourReset}")
endif()

# Enable FXT if StarPU
option(CHAMELEON_ENABLE_TRACING "Enable tracing support" OFF)
if (NOT CHAMELEON_ENABLE_TRACING)
    message("-- ${BoldGreen}CHAMELEON_ENABLE_TRACING is set to OFF, turn it ON to use FxT (with StarPU)${ColourReset}")
endif()
#option(CHAMELEON_USE_EZTRACE "Enable EZTRACE to build modules" OFF)

option(CHAMELEON_RUNTIME_SYNC "Enable synchronous task submission when available to debug the code without parallelism" OFF)
if (CHAMELEON_RUNTIME_SYNC)
    message("-- ${BoldGreen}CHAMELEON_RUNTIME_SYNC is set to ON, turn it OFF to avoid synchronisation in the tasks submission${ColourReset}")
endif()

# Options to enable/disable testings and timings
# ----------------------------------------------
option(CHAMELEON_ENABLE_DOC      "Enable documentation build"  OFF)
if (CHAMELEON_ENABLE_DOC)
    message("-- ${BoldGreen}CHAMELEON_ENABLE_DOC is set to ON, turn it OFF to avoid building docs${ColourReset}")
endif()
option(CHAMELEON_ENABLE_EXAMPLE  "Enable examples build"       ON)
if (CHAMELEON_ENABLE_EXAMPLE)
    message("-- ${BoldGreen}CHAMELEON_ENABLE_EXAMPLE is set to ON, turn it OFF to avoid building examples${ColourReset}")
endif()
option(CHAMELEON_ENABLE_TESTING  "Enable testings build"       ON)
if (CHAMELEON_ENABLE_TESTING)
    message("-- ${BoldGreen}CHAMELEON_ENABLE_TESTING is set to ON, turn it OFF to avoid building testing${ColourReset}")
endif()
option(CHAMELEON_ENABLE_TIMING   "Enable timings build"        ON)
if (CHAMELEON_ENABLE_TIMING)
    message("-- ${BoldGreen}CHAMELEON_ENABLE_TIMING is set to ON, turn it OFF to avoid building timing${ColourReset}")
endif()

# Option to activate or not simulation mode (use Simgrid through StarPU)
# ----------------------------------------------------------------------
if(CHAMELEON_SCHED_STARPU)
    option(CHAMELEON_SIMULATION "Enable simulation mode using Simgrid through StarPU" OFF)
    if (NOT CHAMELEON_SIMULATION)
        message("-- ${BoldGreen}CHAMELEON_SIMULATION is set to OFF, turn it ON to use"
            " SIMULATION mode (only with StarPU compiled with SimGrid)${ColourReset}")
    endif()
    option(CHAMELEON_ENABLE_PRUNING_STATS "Enable pruning statistics" OFF)
    if (NOT CHAMELEON_ENABLE_PRUNING_STATS)
        message("-- ${BoldGreen}CHAMELEON_ENABLE_PRUNING_STATS is set to OFF, turn it ON to build pruning statistics${ColourReset}")
    endif()
    if(CHAMELEON_ENABLE_PRUNING_STATS)
        add_definitions(-DCHAMELEON_ENABLE_PRUNING_STATS)
    endif(CHAMELEON_ENABLE_PRUNING_STATS)
endif()

# Initially we need to generate files for different precisions
# TODO: use this option to avoid generate multiple precisions each time we launch cmake
#option(CHAMELEON_GEN_PREC "Generate source files precisions" ON)
#------------------------------------------------------------------------------

###############################################################################
# Build dependency HQR library #
################################
add_subdirectory(hqr)
include_directories(hqr/include)

###############################################################################
# Look for dependencies #
#########################
set(CHAMELEON_DEP "")

# Check for Thread library
# ------------------------
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
if( THREADS_FOUND AND CMAKE_THREAD_LIBS_INIT )
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
  list(APPEND EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# Add math library to the list of extra
# it normally exists on all common systems provided with a C compiler
if(UNIX OR WIN32)
    find_library(
      M_LIBRARY
      NAMES m
      )
    mark_as_advanced(M_LIBRARY)
    if (M_LIBRARY)
      list(APPEND EXTRA_LIBRARIES "${M_LIBRARY}")
    else()
      message(FATAL_ERROR "Could NOT find libm on your system."
        " Are you sure to a have a C compiler installed?")
    endif()
endif()

# Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
# on Unix systems except Apple ones because it does not exist on it
if(UNIX AND NOT APPLE)
    find_library(
        RT_LIBRARY
        NAMES rt
        )
    mark_as_advanced(RT_LIBRARY)
    if (RT_LIBRARY)
        list(APPEND EXTRA_LIBRARIES "${RT_LIBRARY}")
    else()
        message(FATAL_ERROR "Could NOT find librt on your system")
    endif()
endif()

# If simulation we don't enter in kernel functions so that we don't need to
# link with concerned libraries
if(NOT CHAMELEON_SIMULATION)
    # CHAMELEON depends on CBLAS
    #---------------------------
    find_package(CBLAS REQUIRED)
    if(CBLAS_FOUND)
        if (CBLAS_INCLUDE_DIRS_DEP)
            include_directories(${CBLAS_INCLUDE_DIRS_DEP})
        endif()
        if(CBLAS_CFLAGS_OTHER_DEP)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CBLAS_CFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
        endif()
        if(CBLAS_LDFLAGS_OTHER_DEP)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CBLAS_LDFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        endif()
        if(CBLAS_LIBRARY_DIRS_DEP)
            list(APPEND CMAKE_INSTALL_RPATH "${CBLAS_LIBRARY_DIRS_DEP}")
        endif()
        if (CBLAS_LIBRARIES_DEP)
            list(INSERT CHAMELEON_DEP 0 ${CBLAS_LIBRARIES_DEP})
        endif()
    else()
        if(CHAMELEON_VERBOSE_FIND_PACKAGE)
            if(CBLAS_STANDALONE OR NOT CBLAS_WORKS)
                if (NOT CBLAS_cblas.h_DIRS)
                    Print_Find_Header_Status(cblas cblas.h)
                endif ()
                if (NOT CBLAS_cblas_LIBRARY)
                    Print_Find_Library_Status(cblas libcblas)
                endif ()
            endif()
        else()
            message(WARNING "CBLAS library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        message(FATAL_ERROR "A CBLAS library is required but has not been found")
    endif()

    # CHAMELEON depends on LAPACKE
    #-----------------------------
    find_package(LAPACKE REQUIRED COMPONENTS TMG)
    if(LAPACKE_FOUND)
        if (LAPACKE_INCLUDE_DIRS_DEP)
            include_directories(${LAPACKE_INCLUDE_DIRS_DEP})
        endif()
        if(LAPACKE_CFLAGS_OTHER_DEP)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LAPACKE_CFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
        endif()
        if(LAPACKE_LDFLAGS_OTHER_DEP)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACKE_LDFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        endif()
        if(LAPACKE_LIBRARY_DIRS_DEP)
            list(APPEND CMAKE_INSTALL_RPATH "${LAPACKE_LIBRARY_DIRS_DEP}")
        endif()
        if (LAPACKE_LIBRARIES_DEP)
            list(INSERT CHAMELEON_DEP 0 ${LAPACKE_LIBRARIES_DEP})
        endif()
    else()
        if(CHAMELEON_VERBOSE_FIND_PACKAGE)
            if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
                if (NOT LAPACKE_lapacke.h_DIRS)
                    Print_Find_Header_Status(lapacke lapacke.h)
                endif ()
                if (NOT LAPACKE_lapacke_LIBRARY)
                    Print_Find_Library_Status(lapacke liblapacke)
                endif ()
            endif()
        else()
            message(WARNING "LAPACKE library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        message(FATAL_ERROR "A LAPACKE library is required but has not been found")
    endif()

    # CHAMELEON depends on CUDA/CUBLAS
    #----------------------------------
    if (CHAMELEON_USE_CUDA)

        find_package(CUDA REQUIRED)
        if (CUDA_FOUND)
            message("-- ${Blue}Add definition CHAMELEON_USE_CUDA"
            " - Activate CUDA in Chameleon${ColourReset}")
            libraries_absolute_path(CUDA_LIBRARIES "$ENV{CUDA_ROOT}/lib64")
            set(EXTRA_LIBRARIES_CUDA)
            if(CUDA_VERSION VERSION_LESS "3.0")
                set(CUDA_HOST_COMPILATION_CPP OFF)
            endif(CUDA_VERSION VERSION_LESS "3.0")
            set(CUDA_BUILD_EMULATION OFF)
            if (CUDA_INCLUDE_DIRS)
                include_directories(${CUDA_INCLUDE_DIRS})
            else()
                message(WARNING "CHAMELEON_USE_CUDA requires"
                "\n   CUDA_INCLUDE_DIRS to be found. Be sure you have"
                "\n   cuda headers with your distribution of CUDA.")
            endif()
            if (CUDA_LIBRARIES)
                set(saved_include "${CMAKE_REQUIRED_INCLUDES}")
                set(saved_libs "${CMAKE_REQUIRED_LIBRARIES}")
                set(CMAKE_REQUIRED_INCLUDES  "${CMAKE_REQUIRED_INCLUDES};${CUDA_INCLUDE_DIRS}")
                set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};${CUDA_LIBRARIES}")
                if(CUDA_VERSION VERSION_LESS "4.0")
                    set(CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS 0)
                else()
                    check_function_exists(cuDeviceCanAccessPeer CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS)
                endif()
                set(CMAKE_REQUIRED_INCLUDES  "${saved_include}")
                set(CMAKE_REQUIRED_LIBRARIES "${saved_libs}")
                # Add cublas if found
                if (CUDA_CUBLAS_LIBRARIES)
                    message("-- ${Blue}Add definition CHAMELEON_USE_CUBLAS"
                    " - Use GPU kernels from cuBLAS${ColourReset}")
                    list(APPEND EXTRA_LIBRARIES_CUDA ${CUDA_CUBLAS_LIBRARIES})
                    list(APPEND EXTRA_LIBRARIES_CUDA ${CUDA_LIBRARIES}) # need CUDA libs to link.
                else()
                    message(FATAL_ERROR "CHAMELEON_USE_CUDA requires"
                    "\n   CUDA_CUBLAS_LIBRARIES to be found. Be sure you have"
                    "\n   libcublas with your distribution of CUDA.")
                endif()
            else()
                message(FATAL_ERROR "CHAMELEON_USE_CUDA requires"
                "\n   CUDA_LIBRARIES to be found. Be sure you have"
                "\n   libcuda with your distribution of CUDA.")
            endif()
        endif (CUDA_FOUND)

        list(APPEND EXTRA_LIBRARIES ${EXTRA_LIBRARIES_CUDA})
        foreach (_lib ${EXTRA_LIBRARIES_CUDA})
            if (EXISTS ${_lib})
                get_filename_component(libdir ${_lib} PATH)
                if (EXISTS ${libdir})
                    list(APPEND CMAKE_INSTALL_RPATH "${libdir}")
                endif()
            endif()
        endforeach()

    endif(CHAMELEON_USE_CUDA)

    # CHAMELEON depends on MPI
    #-------------------------
    if (CHAMELEON_USE_MPI)

        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER mpicc)
        endif()
        find_package(MPI REQUIRED)

        if (MPI_C_FOUND)
            message("-- ${Blue}Add definition CHAMELEON_USE_MPI"
            " - Activate MPI in Chameleon${ColourReset}")
            list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
            include_directories( ${MPI_C_INCLUDE_PATH} )
            # tests for intel mpi
            #list(APPEND MPI_C_COMPILE_FLAGS "-mt_mpi")
            #list(APPEND MPI_COMPILE_FLAGS "-mt_mpi")
            if(MPI_C_COMPILE_FLAGS)
                set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
                string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
            endif()
            if(MPI_C_LINK_FLAGS)
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
                string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
            endif()
        endif (MPI_C_FOUND)

        foreach (_lib ${MPI_C_LIBRARIES})
            if (EXISTS ${_lib})
                get_filename_component(libdir ${_lib} PATH)
                if (EXISTS ${libdir})
                    list(APPEND CMAKE_INSTALL_RPATH "${libdir}")
                endif()
            endif()
        endforeach()

    endif (CHAMELEON_USE_MPI)

else (NOT CHAMELEON_SIMULATION)

    message("-- ${BoldBlue}Simulation mode of CHAMELEON is activated (CHAMELEON_SIMULATION=ON)."
    "\n   This mode allows you to simulate execution of algorithms with StarPU compiled with SimGrid."
    "\n   To do so, we provide some perfmodels in the simucore/perfmodels/ directory of CHAMELEON sources."
    "\n   The resulting directory where to find perfmodels is: ${CMAKE_CURRENT_SOURCE_DIR}/simucore/perfmodels."
    "\n   To use these perfmodels, please set your STARPU_HOME environment variable to ${CMAKE_CURRENT_SOURCE_DIR}/simucore/perfmodels."
    "\n   Finally, you need to set your STARPU_HOSTNAME environment variable to the name of the machine to simulate."
    "\n   For example: STARPU_HOSTNAME=mirage."
    "\n   Note that only POTRF kernels with block sizes of 320 or 960 (simple and double precision) on mirage machine are available for now."
    "\n   Database of models is subject to change, it should be enrich in a near future.${ColourReset}")

    message("-- ${Blue}Add definition CHAMELEON_SIMULATION"
    " - Activate simulation mode in Chameleon (to use with StarPU+SimGrid)${ColourReset}")

    add_subdirectory(simucore)

    # Guard against mis-used simulation mode
    if(NOT DEFINED CHAMELEON_USE_CUDA)
        message(FATAL_ERROR "${BoldBlue}"
        "In simulation mode CHAMELEON_USE_CUDA should be set to"
        "   ensure that the user is aware of the version to be used. If Chameleon's"
        "   kernels are available for NVIDIA CUDA GPUs and if the according"
        "   perfmodels are available in STARPU_HOME then use CHAMELEON_USE_CUDA=ON"
        "   else set CHAMELEON_USE_CUDA=OFF."
        "${ColourReset}")
    endif()

    if (NOT CHAMELEON_SCHED_STARPU)
        message(FATAL_ERROR "Simulation mode of Chameleon only works with"
        "\n   the StarPU runtime. Please use CHAMELEON_SCHED_STARPU=ON.")
    endif()

    if (CHAMELEON_ENABLE_EXAMPLE)
        set(CHAMELEON_ENABLE_EXAMPLE OFF)
        message("-- ${BoldBlue}CHAMELEON_ENABLE_EXAMPLE is set to ON but we turn it OFF.")
    endif ()
    if (CHAMELEON_ENABLE_TESTING)
        set(CHAMELEON_ENABLE_TESTING OFF)
        message("-- ${BoldBlue}CHAMELEON_ENABLE_TESTING is set to ON but we turn it OFF."
            "\n   Because we are compiling the simulation mode (CHAMELEON_SIMULATION=ON),"
            "\n   there is no sense in compiling testing drivers that are used to check"
            "\n   numerical correctness of algorithms and kernels.${ColourReset}")
    endif ()

    # Simulation mode: we depend on SimGrid
    find_package(SIMGRID REQUIRED)

    # CHAMELEON depends on MPI
    #-------------------------
    if (CHAMELEON_USE_MPI)

        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpicc")
        endif()
        if(NOT MPI_C_LIBRARIES)
            set(MPI_C_LIBRARIES "${SIMGRID_LIBRARIES}")
        endif()
        if(NOT MPI_C_INCLUDE_PATH)
            set(MPI_C_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
            list(APPEND MPI_C_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
        endif()
        if(NOT MPI_CXX_COMPILER)
            set(MPI_CXX_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpicxx")
        endif()
        if(NOT MPI_CXX_LIBRARIES)
            set(MPI_CXX_LIBRARIES "${SIMGRID_LIBRARIES}")
        endif()
        if(NOT MPI_CXX_INCLUDE_PATH)
            set(MPI_CXX_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
            list(APPEND MPI_CXX_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
        endif()
        if(NOT MPI_Fortran_COMPILER)
            set(MPI_Fortran_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpicxx")
        endif()
        if(NOT MPI_Fortran_LIBRARIES)
            set(MPI_Fortran_LIBRARIES "${SIMGRID_LIBRARIES}")
        endif()
        if(NOT MPI_Fortran_INCLUDE_PATH)
            set(MPI_Fortran_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
            list(APPEND MPI_Fortran_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
        endif()
        find_package(MPI REQUIRED)

        if (MPI_C_FOUND)
            message("-- ${Blue}Add definition CHAMELEON_USE_MPI"
            " - Activate MPI in Chameleon${ColourReset}")
            list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
            include_directories( ${MPI_C_INCLUDE_PATH} )
            if(MPI_C_COMPILE_FLAGS)
                set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
                string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
            endif()
            if(MPI_C_LINK_FLAGS)
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
                string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
            endif()
        endif (MPI_C_FOUND)

    endif (CHAMELEON_USE_MPI)

    foreach (_lib ${SIMGRID_LIBRARIES} ${MPI_C_LIBRARIES})
        if (EXISTS ${_lib})
            get_filename_component(libdir ${_lib} PATH)
            if (EXISTS ${libdir})
                list(APPEND CMAKE_INSTALL_RPATH "${libdir}")
            endif()
        endif()
    endforeach()

endif(NOT CHAMELEON_SIMULATION)


# CHAMELEON depends on a runtime
# -------------------------------
if( CHAMELEON_SCHED_STARPU )

    set(CHAMELEON_STARPU_VERSION "1.3" CACHE STRING "oldest STARPU version desired")

    # create list of components in order to make a single call to find_package(starpu...)
    if(NOT CHAMELEON_SIMULATION)
        set(STARPU_COMPONENT_LIST "HWLOC")
        if(CHAMELEON_USE_CUDA)
            list(APPEND STARPU_COMPONENT_LIST "CUDA")
        endif()
    else()
        set(STARPU_COMPONENT_LIST "SIMGRID")
    endif()
    if(CHAMELEON_USE_MPI)
        list(APPEND STARPU_COMPONENT_LIST "MPI")
    endif()
    if(CHAMELEON_ENABLE_TRACING)
        list(APPEND STARPU_COMPONENT_LIST "FXT")
    endif()

    find_package(STARPU ${CHAMELEON_STARPU_VERSION} REQUIRED
                 COMPONENTS ${STARPU_COMPONENT_LIST})

    # Add definition and include_dir if found
    if ( STARPU_FOUND )
        message("-- ${Blue}Add definition CHAMELEON_SCHED_STARPU"
        " - Activate StarPU in Chameleon${ColourReset}")
        if (STARPU_INCLUDE_DIRS_DEP)
            include_directories(${STARPU_INCLUDE_DIRS_DEP})
            set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS_DEP}")
        endif()
        if(STARPU_CFLAGS_OTHER_DEP)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${STARPU_CFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
            set(CMAKE_REQUIRED_FLAGS "${STARPU_CFLAGS_OTHER_DEP}")
        endif()
        if(STARPU_LDFLAGS_OTHER_DEP)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STARPU_LDFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
            set(CMAKE_REQUIRED_LDFLAGS "${STARPU_LDFLAGS_OTHER_DEP}")
        endif()
        if(STARPU_LIBRARY_DIRS_DEP)
            list(APPEND CMAKE_INSTALL_RPATH "${STARPU_LIBRARY_DIRS_DEP}")
        endif()
        if (STARPU_LIBRARIES_DEP)
            list(INSERT CHAMELEON_DEP 0 ${STARPU_LIBRARIES_DEP})
            set(CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES_DEP}")
        endif()
        string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
        check_function_exists(starpu_data_idle_prefetch_on_node HAVE_STARPU_IDLE_PREFETCH)
        if ( HAVE_STARPU_IDLE_PREFETCH )
            message("-- ${Blue}Add definition HAVE_STARPU_IDLE_PREFETCH${ColourReset}")
        endif()
        check_function_exists(starpu_iteration_push HAVE_STARPU_ITERATION_PUSH)
        if ( HAVE_STARPU_ITERATION_PUSH )
            message("-- ${Blue}Add definition HAVE_STARPU_ITERATION_PUSH${ColourReset}")
        endif()
        check_function_exists(starpu_data_wont_use HAVE_STARPU_DATA_WONT_USE)
        if ( HAVE_STARPU_DATA_WONT_USE )
            message("-- ${Blue}Add definition HAVE_STARPU_DATA_WONT_USE${ColourReset}")
        endif()
        check_function_exists(starpu_data_set_coordinates HAVE_STARPU_DATA_SET_COORDINATES)
        if ( HAVE_STARPU_DATA_SET_COORDINATES )
            message("-- ${Blue}Add definition HAVE_STARPU_DATA_SET_COORDINATES${ColourReset}")
        endif()
        check_function_exists(starpu_malloc_on_node_set_default_flags HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS)
        if ( HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS )
            message("-- ${Blue}Add definition HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS${ColourReset}")
        endif()
        check_function_exists(starpu_data_set_ooc_flag HAVE_STARPU_DATA_SET_OOC_FLAG)
        if ( HAVE_STARPU_DATA_SET_OOC_FLAG )
            message("-- ${Blue}Add definition HAVE_STARPU_DATA_SET_OOC_FLAG${ColourReset}")
        endif()
        check_function_exists(starpu_mpi_init_conf HAVE_STARPU_MPI_INIT_CONF)
        if ( HAVE_STARPU_MPI_INIT_CONF )
            message("-- ${Blue}Add definition HAVE_STARPU_MPI_INIT_CONF${ColourReset}")
        endif()
        check_function_exists(starpu_mpi_data_migrate HAVE_STARPU_MPI_DATA_MIGRATE)
        if ( HAVE_STARPU_MPI_DATA_MIGRATE )
            message("-- ${Blue}Add definition HAVE_STARPU_MPI_DATA_MIGRATE${ColourReset}")
        elseif(CHAMELEON_USE_MIGRATE)
            set(CHAMELEON_USE_MIGRATE "OFF")
            message("-- ${Blue}CHAMELEON_USE_MIGRATE is turned OFF because starpu_mpi_data_migrate not found${ColourReset}")
        endif()
        if(CHAMELEON_ENABLE_TRACING)
            # check if fxt profiling is accessible in starpu and activate it in chameleon
            check_function_exists(starpu_fxt_start_profiling HAVE_STARPU_FXT_PROFILING)
            if ( HAVE_STARPU_FXT_PROFILING )
                message("-- ${Blue}Add definition HAVE_STARPU_FXT_PROFILING"
                " - Activate FxT profiling through StarPU${ColourReset}")
            else()
                message("-- ${Red}Looking for starpu with fxt"
                " - starpu_fxt_start_profiling() test fails in StarPU${ColourReset}")
                message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
            endif()
        endif()
        if (CHAMELEON_USE_MPI)
            # Check if a specific function exist
            check_function_exists(starpu_mpi_data_register_comm HAVE_STARPU_MPI_DATA_REGISTER)
            if ( HAVE_STARPU_MPI_DATA_REGISTER )
                message("-- ${Blue}Add definition HAVE_STARPU_MPI_DATA_REGISTER - Activate"
                " use of starpu_mpi_data_register() in Chameleon with StarPU${ColourReset}")
            endif()
            check_function_exists(starpu_mpi_comm_rank HAVE_STARPU_MPI_COMM_RANK)
            if ( HAVE_STARPU_MPI_COMM_RANK )
                message("-- ${Blue}Add definition HAVE_STARPU_MPI_COMM_RANK - Activate"
                " use of starpu_mpi_comm_rank() in Chameleon with StarPU${ColourReset}")
            endif()
            check_function_exists(starpu_mpi_cached_receive HAVE_STARPU_MPI_CACHED_RECEIVE)
            if ( HAVE_STARPU_MPI_CACHED_RECEIVE )
                message("-- ${Blue}Add definition HAVE_STARPU_MPI_CACHED_RECEIVE${ColourReset}")
            endif()
            check_function_exists(starpu_mpi_comm_get_attr HAVE_STARPU_MPI_COMM_GET_ATTR)
            if ( HAVE_STARPU_MPI_COMM_GET_ATTR )
                message("-- ${Blue}Add definition HAVE_STARPU_MPI_COMM_GET_ATTR${ColourReset}")
            endif()
        endif()
    else ( STARPU_FOUND )
        if(CHAMELEON_VERBOSE_FIND_PACKAGE)
            if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                if (NOT HWLOC_hwloc.h_DIRS)
                    Print_Find_Header_Status(hwloc hwloc.h)
                endif ()
                if (NOT HWLOC_hwloc_LIBRARY)
                    Print_Find_Library_Status(hwloc libhwloc)
                endif ()
            endif()
            if(CHAMELEON_ENABLE_TRACING AND (NOT FXT_FOUND OR NOT FXT_LIBRARIES))
                if (NOT FXT_fxt.h_DIRS)
                    Print_Find_Header_Status(fxt fxt.h)
                endif ()
                if (NOT FXT_fxt_LIBRARY)
                    Print_Find_Library_Status(fxt libfxt)
                endif ()
            endif()
            if(CHAMELEON_SIMULATION AND (NOT SIMGRID_FOUND OR NOT SIMGRID_LIBRARIES))
                if (NOT SIMGRID_simgrid.h_DIRS)
                    Print_Find_Header_Status(simgrid simgrid.h)
                endif ()
                if (NOT SIMGRID_simgrid_LIBRARY)
                    Print_Find_Library_Status(simgrid libsimgrid)
                endif ()
            endif()
            if( (NOT STARPU_SHM_FOUND) OR (NOT STARPU_SHM_LIBRARIES) OR
                ( STARPU_LOOK_FOR_MPI AND (NOT STARPU_MPI_FOUND OR NOT STARPU_MPI_LIBRARIES) )
              )
                foreach(starpu_hdr ${STARPU_hdrs_to_find})
                    if (NOT STARPU_${starpu_hdr}_INCLUDE_DIRS)
                        Print_Find_Header_Status(starpu ${starpu_hdr})
                    endif ()
                endforeach()
                if(STARPU_VERSION_STRING)
                    foreach(starpu_lib ${STARPU_libs_to_find})
                        if (NOT STARPU_${starpu_lib}_LIBRARY)
                            Print_Find_Library_Status(starpu ${starpu_lib})
                        endif ()
                    endforeach()
                endif ()
            endif ()
        else(CHAMELEON_VERBOSE_FIND_PACKAGE)
            message(WARNING "StarPU library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif(CHAMELEON_VERBOSE_FIND_PACKAGE)
        if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
            message(FATAL_ERROR "hwloc library is required but has not been found")
        endif()
        if(CHAMELEON_SIMULATION AND (NOT SIMGRID_FOUND OR NOT SIMGRID_LIBRARIES))
            message(FATAL_ERROR "SimGrid library is required but has not been found")
        endif()
        if(CHAMELEON_ENABLE_TRACING AND (NOT FXT_FOUND OR NOT FXT_LIBRARIES))
            message(FATAL_ERROR "FxT library is required but has not been found")
        endif()
        if( (NOT STARPU_SHM_FOUND) OR (NOT STARPU_SHM_LIBRARIES) OR
            ( STARPU_LOOK_FOR_MPI AND (NOT STARPU_MPI_FOUND OR NOT STARPU_MPI_LIBRARIES) )
          )
            message(FATAL_ERROR "StarPU library is required but has not been found")
        endif()
    endif ( STARPU_FOUND )

endif( CHAMELEON_SCHED_STARPU )

if( CHAMELEON_SCHED_PARSEC )

    # create list of components in order to make a single call to find_package(starpu...)
    set(PARSEC_COMPONENT_LIST "HWLOC")
    if(CHAMELEON_USE_MPI)
        list(APPEND PARSEC_COMPONENT_LIST "MPI")
    endif()
    if(CHAMELEON_USE_CUDA)
        list(APPEND PARSEC_COMPONENT_LIST "CUDA")
    endif()
    # TODO: Add a CHAMELEON_WITH_PROFILING option that enables Fxt for StarPU, or PAPI for PaRSEC
    #if(CHAMELEON_WITH_PROFILING)
    #    list(APPEND PARSEC_COMPONENT_LIST "PAPI")
    #endif()
    find_package(PARSEC COMPONENTS ${PARSEC_COMPONENT_LIST})
    if(PARSEC_FOUND)
        message("-- ${Blue}Add definition CHAMELEON_SCHED_PARSEC"
        " - Activate PaRSEC in Chameleon${ColourReset}")
      if (PARSEC_INCLUDE_DIRS_DEP)
          include_directories(${PARSEC_INCLUDE_DIRS_DEP})
      endif()
      if(PARSEC_CFLAGS_OTHER_DEP)
          set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PARSEC_CFLAGS_OTHER_DEP}")
          string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
      endif()
      if(PARSEC_LDFLAGS_OTHER_DEP)
          set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PARSEC_LDFLAGS_OTHER_DEP}")
          string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
      endif()
      if(PARSEC_LIBRARY_DIRS_DEP)
          list(APPEND CMAKE_INSTALL_RPATH "${PARSEC_LIBRARY_DIRS_DEP}")
      endif()
      if (PARSEC_LIBRARIES_DEP)
          list(INSERT CHAMELEON_DEP 0 ${PARSEC_LIBRARIES_DEP})
      endif()
      string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
      # to test if symbols exist, e.g.
      # set(CMAKE_REQUIRED_INCLUDES "${PARSEC_INCLUDE_DIRS_DEP}")
      # set(CMAKE_REQUIRED_FLAGS "${PARSEC_CFLAGS_OTHER_DEP}")
      # set(CMAKE_REQUIRED_LDFLAGS "${PARSEC_LDFLAGS_OTHER_DEP}")
      # set(CMAKE_REQUIRED_LIBRARIES "${PARSEC_LIBRARIES_DEP}")
      # check_function_exists(parsec_init HAVE_PARSEC_INIT)
    else ( PARSEC_FOUND )
        if(CHAMELEON_VERBOSE_FIND_PACKAGE)
            if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                if (NOT HWLOC_hwloc.h_DIRS)
                    Print_Find_Header_Status(hwloc hwloc.h)
                endif ()
                if (NOT HWLOC_hwloc_LIBRARY)
                    Print_Find_Library_Status(hwloc libhwloc)
                endif ()
            endif()
            if( (NOT PARSEC_FOUND) OR (NOT PARSEC_LIBRARIES_DEP) )
              foreach(parsec_hdr ${PARSEC_hdrs_to_find})
                    if (NOT PARSEC_${parsec_hdr}_INCLUDE_DIRS)
                        Print_Find_Header_Status(parsec ${parsec_hdr})
                    endif ()
                endforeach()
                if(PARSEC_VERSION_STRING)
                    foreach(parsec_lib ${PARSEC_libs_to_find})
                        if (NOT PARSEC_${parsec_lib}_LIBRARY)
                            Print_Find_Library_Status(parsec ${parsec_lib})
                        endif ()
                    endforeach()
                endif()
            endif()
        else(CHAMELEON_VERBOSE_FIND_PACKAGE)
            message(WARNING "PaRSEC library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif(CHAMELEON_VERBOSE_FIND_PACKAGE)
        if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
            message(FATAL_ERROR "hwloc library is required but has not been found")
        endif()
        if( (NOT PARSEC_FOUND) OR (NOT PARSEC_LIBRARIES_DEP) )
            message(FATAL_ERROR "PaRSEC library is required but has not been found")
        endif()
    endif ( PARSEC_FOUND )

endif( CHAMELEON_SCHED_PARSEC )

      # try to find quark runtime
if( CHAMELEON_SCHED_QUARK )

    find_package(QUARK COMPONENTS HWLOC)

    # Add definition and include_dir if found
    if(QUARK_FOUND)
        message("-- ${Blue}Add definition CHAMELEON_SCHED_QUARK"
        " - Activate QUARK in Chameleon${ColourReset}")
        if (QUARK_INCLUDE_DIRS_DEP)
            include_directories(${QUARK_INCLUDE_DIRS_DEP})
        endif()
        if(QUARK_CFLAGS_OTHER_DEP)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${QUARK_CFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
        endif()
        if(QUARK_LDFLAGS_OTHER_DEP)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${QUARK_LDFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        endif()
        if(QUARK_LIBRARY_DIRS_DEP)
            list(APPEND CMAKE_INSTALL_RPATH "${QUARK_LIBRARY_DIRS_DEP}")
        endif()
        if (QUARK_LIBRARIES_DEP)
            list(INSERT CHAMELEON_DEP 0 ${QUARK_LIBRARIES_DEP})
        endif()
    else(QUARK_FOUND)
        if(CHAMELEON_VERBOSE_FIND_PACKAGE)
            if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                if (NOT HWLOC_hwloc.h_DIRS)
                    Print_Find_Header_Status(hwloc hwloc.h)
                endif ()
                if (NOT HWLOC_hwloc_LIBRARY)
                    Print_Find_Library_Status(hwloc libhwloc)
                endif ()
            endif()
            if (NOT QUARK_quark.h_DIRS)
                Print_Find_Header_Status(quark quark.h)
            endif ()
            if (NOT QUARK_quark_LIBRARY)
                Print_Find_Library_Status(quark libquark)
            endif ()
        else()
            message(WARNING "QUARK library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
            message(FATAL_ERROR "hwloc library is required but has not been found")
        endif()
        message(FATAL_ERROR "QUARK library is required but has not been found")
    endif(QUARK_FOUND)

endif()

if( CHAMELEON_SCHED_OPENMP )
    find_package(OpenMP REQUIRED)

    if ( OPENMP_FOUND )
        message("-- ${Blue}Add definition CHAMELEON_SCHED_OPENMP"
        " - Activate OpenMP in Chameleon${ColourReset}")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
        # TODO: check for compiler support for "_Pragma"
    else ( OPENMP_FOUND )
      message(FATAL_ERROR "Something went wrong when finding OpenMP")
    endif ( OPENMP_FOUND )
endif( CHAMELEON_SCHED_OPENMP )

# getopt
unset(CMAKE_REQUIRED_LIBRARIES)
check_include_files(getopt.h CHAMELEON_HAVE_GETOPT_H)
if (CHAMELEON_HAVE_GETOPT_H)
  check_function_exists(getopt_long CHAMELEON_HAVE_GETOPT_LONG)
endif()

# Add option to exploit cublas API v2
# -----------------------------------
cmake_dependent_option(CHAMELEON_USE_CUBLAS_V2
  "Enable cublas API v2" ON
  "CHAMELEON_USE_CUDA;CHAMELEON_SCHED_STARPU" OFF)

list(REMOVE_DUPLICATES CMAKE_C_FLAGS)
string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
list(REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS)
string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
# Fix a problem on Mac OS X when building shared libraries
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_SHARED_LINKER_FLAGS "-undefined dynamic_lookup")
endif()

# Save extra dependencies (all required links)
list(APPEND CHAMELEON_DEP ${EXTRA_LIBRARIES})
# We remove duplicate dependencies ensuring that redundant libraries are placed last
list(REVERSE CHAMELEON_DEP)
list(REMOVE_DUPLICATES CHAMELEON_DEP)
list(REVERSE CHAMELEON_DEP)

# we add chameleon's own library
set(CHAMELEON_LIBRARIES "chameleon")
if (NOT CHAMELEON_SIMULATION)
  list(APPEND CHAMELEON_LIBRARIES coreblas)
  if (CHAMELEON_USE_CUDA)
    list(APPEND CHAMELEON_LIBRARIES cudablas)
  endif()
endif()
if (CHAMELEON_SCHED_QUARK)
  list(APPEND CHAMELEON_LIBRARIES chameleon_quark)
endif()
if (CHAMELEON_SCHED_PARSEC)
  list(APPEND CHAMELEON_LIBRARIES chameleon_parsec)
endif()
if (CHAMELEON_SCHED_STARPU)
  list(APPEND CHAMELEON_LIBRARIES chameleon_starpu)
endif()
if (CHAMELEON_SCHED_OPENMP)
  list(APPEND CHAMELEON_LIBRARIES chameleon_openmp)
endif()
list(APPEND CHAMELEON_LIBRARIES hqr)
set(CHAMELEON_LIBRARIES_DEP ${CHAMELEON_LIBRARIES} ${CHAMELEON_DEP})

# Check for the subdirectories
# ----------------------------

set(CHAMELEON_SOURCES_TARGETS "" CACHE INTERNAL "List of targets of sources")

# include headers
# ---------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

if (NOT CHAMELEON_SIMULATION)
  ###############################################################################
  # Coreblas library (kernels for CPUs) #
  #######################################
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/coreblas/include)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/coreblas/include)
  add_subdirectory(coreblas)

  ###############################################################################
  # Cudablas library (kernels for CUDAs) #
  ########################################
  if(CHAMELEON_USE_CUDA)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cudablas/include)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/cudablas/include)
    add_subdirectory(cudablas)
  endif()

endif()

###############################################################################
# Main library #
################
add_subdirectory(include)
add_subdirectory(control)
add_subdirectory(compute)

###############################################################################
# Runtimes libraries #
######################
add_subdirectory(runtime)

###############################################################################
# Executables and Tests #
#########################
enable_testing()
include(CTest)

# Necessary to compile C executables linked with ifort
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    message(STATUS "Add -nofor_main to the Fortran linker (Intel compiler)")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -nofor_main")
endif()

# Testing executables
if(CHAMELEON_ENABLE_EXAMPLE AND NOT CHAMELEON_SIMULATION)
  add_subdirectory(example)
endif(CHAMELEON_ENABLE_EXAMPLE AND NOT CHAMELEON_SIMULATION)

# Testing executables
if(CHAMELEON_ENABLE_TESTING AND NOT CHAMELEON_SIMULATION)
  add_subdirectory(testing)
  add_subdirectory(new-testing)
endif(CHAMELEON_ENABLE_TESTING AND NOT CHAMELEON_SIMULATION)

# Timing executables
if(CHAMELEON_ENABLE_TIMING)
  add_subdirectory(timing)
endif(CHAMELEON_ENABLE_TIMING)
#------------------------------------------------------------------------------

# Define a target which gathers all targets of sources
# useful for generating the doc without compiling the binaries (objects, libraries)
add_custom_target(chameleon_all_sources ALL DEPENDS ${CHAMELEON_SOURCES_TARGETS})

###############################################################################
# Documentation #
#################

if(CHAMELEON_ENABLE_DOC)
    add_subdirectory(doc)
    if (DOXYGEN_EXECUTABLE AND EMACS_COMPILER AND TEX_COMPILER)
        add_custom_target(doc ALL DEPENDS doxygen-out doc-html-users_guide doc-pdf-users_guide)
    endif()
endif()
#------------------------------------------------------------------------------


###############################################################################
# Make cmake include starpu_simgrid_wrap.h (which #define main and thus makes #
# link tests fail) only after tests.                                          #
###############################################################################
if( CHAMELEON_SCHED_STARPU )
    if ( STARPU_FOUND )
        if (CHAMELEON_SIMULATION)
            list(APPEND CMAKE_REQUIRED_FLAGS "-include" "starpu_simgrid_wrap.h")
        endif()
    endif()
endif()

###############################################################################
# Config files (pkg_config) #
#############################

# Create .pc file
# ---------------
generate_pkgconfig_file()
#------------------------------------------------------------------------------

###############################################################################
# Release source #
##################
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_NAME "chameleon")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.org)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Chameleon is a C library providing parallel algorithms to perform BLAS/LAPACK "
                                      " operations exploiting fully modern architectures.")
set(CHAMELEON_VERSION_PATCH ${CHAMELEON_VERSION_MICRO})
set(CPACK_PACKAGE_VERSION "${CHAMELEON_VERSION_MAJOR}.${CHAMELEON_VERSION_MINOR}.${CHAMELEON_VERSION_PATCH}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
set(CPACK_PACKAGE_VENDOR "Inria")
set(CPACK_PACKAGE_CONTACT "morse-devel@lists.gforge.inria.fr")
set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.org)
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENCE.txt)
set(CPACK_SOURCE_IGNORE_FILES "/plasma-conversion/;/build*;.cproject;.settings;.dir-locals.el;.project;.pydevproject;/.git*")
include(CPack)

###############################################################################
# Print Options #
#################
include(PrintOpts)

if (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files")
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} which is a default system path."
    " You may want to change it: set the CMAKE_INSTALL_PREFIX variable to do so${ColourReset}")
else()
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}${ColourReset}")
endif()

###
### END CMakeLists.txt
###
