# --------------------------------------------------------------------------- #
#    CMake file for BinaryKnapsackBlock/tools                                 #
# --------------------------------------------------------------------------- #

# This is needed for setting the runtime path when installing.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Since we are using the block factory, objects from linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.
if (BUILD_SHARED_LIBS)
    if (UNIX AND (NOT APPLE))
        add_link_options("-Wl,--no-as-needed")
    endif ()
else ()
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
