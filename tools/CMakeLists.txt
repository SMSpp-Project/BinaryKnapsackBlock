# --------------------------------------------------------------------------- #
#    CMake file for BinaryKnapsackBlock/tools                                 #
# --------------------------------------------------------------------------- #

# This is needed for setting the runtime path when installing.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The classes of the modules register themselves in the factories from
# static objects, which nothing here names: a shared library is kept by
# what its module asks of the linker [see the CMakeLists.txt of the
# module], while an archive has to be taken whole.
if (NOT BUILD_SHARED_LIBS)
    if (MSVC)
        add_link_options("/WHOLEARCHIVE" "/FORCE:MULTIPLE")
    else () # Unix
        if (APPLE)
            add_link_options("-Wl,-all_load")
        else ()
            add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
        endif ()
    endif ()
endif ()

# ----- bkbench ------------------------------------------------------------ #
# in-process benchmark over the large public Pisinger/Jooken archives (the
# default solver is the efficient CoreDPBinaryKnapsackSolver)
add_executable(bkbench bkbench.cpp)
target_link_libraries(bkbench PRIVATE SMS++::BinaryKnapsackBlock)

# ----- bk2nc4 -------------------------------------------------------------- #
# the textual instances of the public archives, in the Pisinger/Jooken format
# or in the native one, written as netCDF files, which is what whoever reads a
# Block rather than a knapsack needs
add_executable(bk2nc4 bk2nc4.cpp)
target_link_libraries(bk2nc4 PRIVATE SMS++::BinaryKnapsackBlock)

# ----- Install instructions ------------------------------------------------ #
include(GNUInstallDirs)
install(TARGETS bkbench bk2nc4
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# --------------------------------------------------------------------------- #
