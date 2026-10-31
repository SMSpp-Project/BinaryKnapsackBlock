# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ BinaryKnapsackBlock                            #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.21)
cmake_policy(VERSION 3.15)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeriveVersion.cmake)
smspp_derive_version(SMSPP_MODULE_VERSION)

project(BinaryKnapsackBlock
        VERSION ${SMSPP_MODULE_VERSION}
        DESCRIPTION "SMS++ BinaryKnapsackBlock module"
        HOMEPAGE_URL https://gitlab.com/smspp/binaryknapsackblock
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})

# C-safe form of modName for the preprocessor macros in <modName>Config.h
string(MAKE_C_IDENTIFIER "${modName}" modMacro)
set(modNamespace "SMS++")

# Find out if it's being called by the umbrella.
get_directory_property(hasParent PARENT_DIRECTORY)

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
# This creates an ENABLE_TESTING option (default: ON) and enables the testing
# using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
include(CTest)

# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

# ----- Requirements -------------------------------------------------------- #
# If it's not being called by the umbrella, we need to
# look for the system-installed libraries.
if (NOT hasParent)
    find_package(SMS++ REQUIRED)
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- RPATH --------------------------------------------------------------- #
# The installed libraries and executables look for the libraries of the same
# installation relative to their own location, so that the installed tree can
# be moved as a whole; the directories of the external libraries are kept.
include(GNUInstallDirs)
if (NOT CMAKE_INSTALL_RPATH)
    if (APPLE)
        set(rpathOrigin "@loader_path")
    else ()
        set(rpathOrigin "$ORIGIN")
    endif ()
    file(RELATIVE_PATH binToLib ${CMAKE_INSTALL_FULL_BINDIR}
            ${CMAKE_INSTALL_FULL_LIBDIR})
    set(CMAKE_INSTALL_RPATH "${rpathOrigin}" "${rpathOrigin}/${binToLib}")
endif ()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_20)

# The name of the shared library carries the version of its binary interface,
# which any minor release may change while the major version is 0, and then
# any major release.
if (PROJECT_VERSION_MAJOR EQUAL 0)
    set(soVersion ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
else ()
    set(soVersion ${PROJECT_VERSION_MAJOR})
endif ()
set_target_properties(${modName} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${soVersion})

# The classes of the module register themselves in the factory from a static
# initialiser, and a linker that drops what looks unused takes the
# registration away with it: the target tells whoever links it to keep the
# symbol that forces the module in, or the library as a whole where naming
# the symbol is not enough.
if (BUILD_SHARED_LIBS)
    if (MSVC)
        target_link_options(${modName} INTERFACE
                "/INCLUDE:SMSpp_force_load_${modName}")
    elseif (UNIX AND (NOT APPLE))
        # the ELF linker decides by the symbols a regular object asks for, so
        # naming the symbol is not enough and the library has to be kept
        target_link_options(${modName} INTERFACE "LINKER:--no-as-needed")
    elseif (APPLE)
        # ld64 keeps a dylib that provides a symbol to the link, one asked
        # for by -u included, even under -dead_strip_dylibs
        target_link_options(${modName} INTERFACE
                "LINKER:-u,_SMSpp_force_load_${modName}")
    endif ()
endif ()

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
# Note: do not GLOB files here.
target_sources(${modName} PRIVATE
        src/BinaryKnapsackBlock.cpp
        src/BinaryKnapsackSolver.cpp
        src/DPBinaryKnapsackSolver.cpp
        src/CoreDPBinaryKnapsackSolver.cpp
        src/GreedyRelaxationBinaryKnapsackSolver.cpp)

# ----- RECORD (optional) --------------------------------------------------- #
# RECORDBinaryKnapsackSolver hands the integer core to RECORD
# (gitlab.com/renanfernandofranco/record), which is distributed as a single
# source file with a main(): given a checkout in RECORD_ROOT, the file is
# copied at configure time into record-shim/RECORD.hpp with the main() cut
# away (the main() reads a bounded knapsack from the standard input, while
# the bridge sets up the 0-1 problem itself), and compiled alone in
# src/RECORDBridge.cpp, the only file that sees its global names. Without a
# RECORD_ROOT, the checkout is looked for in the environment variable of the
# same name and where INSTALL.sh puts it; RECORD uses POSIX clocks, which
# MSVC does not have.
find_path(RECORD_ROOT NAMES source/RECORD.cpp
        HINTS ENV RECORD_ROOT PATHS /opt/record /Library/record
        DOC "Checkout of RECORD (optional)" NO_DEFAULT_PATH)
if (RECORD_ROOT AND EXISTS ${RECORD_ROOT}/source/RECORD.cpp AND NOT MSVC)
    file(READ ${RECORD_ROOT}/source/RECORD.cpp _recordSrc)
    string(FIND "${_recordSrc}" "int main (" _recordMain)
    if (_recordMain GREATER -1)
        string(SUBSTRING "${_recordSrc}" 0 ${_recordMain} _recordSrc)
    endif ()
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/record-shim/RECORD.hpp
            "${_recordSrc}")
    target_sources(${modName} PRIVATE
            src/RECORDBridge.cpp
            src/RECORDBinaryKnapsackSolver.cpp)
    set_source_files_properties(src/RECORDBridge.cpp PROPERTIES
            INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/record-shim)
    if (NOT MSVC)               # its simd pragmas
        set_source_files_properties(src/RECORDBridge.cpp PROPERTIES
                COMPILE_OPTIONS "-fopenmp-simd")
    endif ()
    message(STATUS "${modName}: RECORD found in ${RECORD_ROOT}.")
endif ()

# ----- COMBO (optional) ---------------------------------------------------- #
# COMBOBinaryKnapsackSolver hands the integer core to COMBO, which is
# distributed by its authors as combo.c and combo.h, for academic or
# non-commercial use only: given a directory holding them in COMBO_ROOT,
# combo.c is copied at configure time into combo-shim/ with its product type
# widened from double to a 128-bit integer (exact determinants of profits
# and weights), and compiled alone in src/COMBOGuard.c, which turns the
# exit() of COMBO on an internal error into an error code. Without a
# COMBO_ROOT, the directory is looked for in the environment variable of the
# same name and in /opt/combo or /Library/combo. The 128-bit integer is not
# available with MSVC.
find_path(COMBO_ROOT NAMES combo.c
        HINTS ENV COMBO_ROOT PATHS /opt/combo /Library/combo
        DOC "Directory with combo.c and combo.h (optional)" NO_DEFAULT_PATH)
if (COMBO_ROOT AND EXISTS ${COMBO_ROOT}/combo.c AND
        EXISTS ${COMBO_ROOT}/combo.h AND NOT MSVC)
    file(READ ${COMBO_ROOT}/combo.c _comboSrc)
    string(REGEX REPLACE "typedef[ \t]+double[ \t]+prod;"
            "typedef __int128 prod;" _comboWide "${_comboSrc}")
    if (_comboWide STREQUAL _comboSrc)
        message(WARNING "${modName}: the type prod of ${COMBO_ROOT}/combo.c "
                "is not the expected double, COMBO is not built.")
    else ()
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/combo-shim/combo.c
                "${_comboWide}")
        configure_file(${COMBO_ROOT}/combo.h
                ${CMAKE_CURRENT_BINARY_DIR}/combo-shim/combo.h COPYONLY)
        target_sources(${modName} PRIVATE
                src/COMBOGuard.c
                src/COMBOBridge.cpp
                src/COMBOBinaryKnapsackSolver.cpp)
        set_source_files_properties(src/COMBOGuard.c src/COMBOBridge.cpp
                PROPERTIES INCLUDE_DIRECTORIES
                ${CMAKE_CURRENT_BINARY_DIR}/combo-shim)
        set_source_files_properties(src/COMBOGuard.c PROPERTIES
                COMPILE_OPTIONS "-w")    # its K&R-era code
        message(STATUS "${modName}: COMBO found in ${COMBO_ROOT}.")
    endif ()
endif ()

# When using target_include_directories(), PUBLIC means that any targets
# that link to this target also need that include directory.
# Other options are PRIVATE (only affect the current target, not dependencies),
# and INTERFACE (only needed for dependencies).

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake).
# This means that if a target finds this library in its build directory
# will look into the BUILD_INTERFACE path for its headers, if it finds it
# installed in the system will look into the INSTALL_INTERFACE path.
target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
target_link_libraries(${modName} PUBLIC ${modNamespace}::SMS++)

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- C++20 module -------------------------------------------------------- #
# The module smspp.binaryknapsackblock (src/BinaryKnapsackBlock.cppm) exports
# the names of the headers and re-exports smspp; like the module smspp of the
# core, it is built when SMSpp_BUILD_MODULE is ON, which needs CMake 3.28 or
# later, the Ninja or Visual Studio generator and a compiler that scans the
# dependencies of the modules.
# Where these are missing, or a module it imports is not there, the option
# is turned off with a warning [see CheckModule.cmake of SMS++].
smspp_check_module(${modNamespace}::SMS++_module)
if (SMSpp_BUILD_MODULE)
    # the scan of the module dependencies, and $<LINK_ONLY> honoured also for
    # the target itself, which cmake_policy(VERSION 3.15) above leaves off
    cmake_policy(SET CMP0155 NEW)
    cmake_policy(SET CMP0131 NEW)

    add_library(${modName}_module)
    target_sources(${modName}_module PUBLIC
            FILE_SET CXX_MODULES
            BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
            FILES src/${modName}.cppm)
    target_compile_features(${modName}_module PUBLIC cxx_std_20)
    # the interface takes the include paths and the definitions of the
    # library, but not its compile options (e.g., the -fopenmp of a source,
    # which no header needs), so that it agrees with the module it imports
    target_include_directories(${modName}_module PRIVATE
            $<TARGET_PROPERTY:${modName},INTERFACE_INCLUDE_DIRECTORIES>)
    target_compile_definitions(${modName}_module PRIVATE
            $<TARGET_PROPERTY:${modName},INTERFACE_COMPILE_DEFINITIONS>)
    target_link_libraries(${modName}_module PUBLIC $<LINK_ONLY:${modName}>
            ${modNamespace}::SMS++_module)
    add_library(${modNamespace}::${modName}_module ALIAS ${modName}_module)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
add_subdirectory(tools)

if (BUILD_TESTING)
    add_subdirectory(test)
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modName}Targets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # .dll
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # .so .dylib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}) # .lib

# Install the module, whose interface is compiled again by each importer
if (SMSpp_BUILD_MODULE)
    install(TARGETS ${modName}_module
            EXPORT ${modName}Targets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            FILE_SET CXX_MODULES
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}/modules)
endif ()

# Install the headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the README and LICENSE files.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/README.md
        ${CMAKE_CURRENT_LIST_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${modName})

# ----- Add the build tree to RPATH------------------------------------------ #
# Add the binary directory to RPATH, so executables that use
# dynamic loading will look into it.
if (hasParent)
    set(SMSPP_BUILD_RPATH
            "${SMSPP_BUILD_RPATH};${CMAKE_CURRENT_BINARY_DIR}"
            PARENT_SCOPE)
endif ()

# ----- Packaging support --------------------------------------------------- #
# The following commands allow one to build a proper package using CPack.
# See: https://cmake.org/cmake/help/latest/manual/cpack.1.html

# Using this with the umbrella project is too messy.
if (NOT hasParent)
    set(CPACK_PACKAGE_VENDOR "SMS++ Team")
    set(CPACK_PACKAGE_CONTACT "Donato Meoli")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_LIST_DIR}/README.md")
    set(CPACK_SOURCE_IGNORE_FILES
            \\\\.git/
            \\\\.gitignore$
            \\\\.gitattributes$)
    include(CPack)
endif ()

# ----- Data ---------------------------------------------------------------- #
# The instances the tests read are, for each <fmt> of DATA_FORMATS, the
# archive <fmt>.tgz in the GitLab Package Registry of the project
# DATA_PROJECT, extracted into data/<fmt>. Nothing of this is built with the
# module: extract_<DATA_PREFIX>_<fmt> downloads the archive (by building
# download_<DATA_PREFIX>_<fmt>) and extracts it, and it is built when the
# tests that read the instances are run, by a test that sets up the fixture
# <DATA_PREFIX>_<fmt>_extracted they require.
#
# A published version is never replaced: new data go out under a new
# version, written here first [see the upload scripts in data/], and the
# archive and the marker of its extraction carry it in their name, so that a
# tree holding the data of an older version (a cache of the CI, or a clone
# that has been extracted before) downloads and extracts again instead of
# taking them for the current ones.
set(DATA_VERSION 2026-06-24)
set(DATA_PROJECT 24425718)
set(DATA_PREFIX bk)
set(DATA_FORMATS txt)

foreach (fmt ${DATA_FORMATS})
    set(url "https://gitlab.com/api/v4/projects/${DATA_PROJECT}/packages/\
generic/${fmt}/${DATA_VERSION}/${fmt}.tgz")
    set(archive ${CMAKE_CURRENT_SOURCE_DIR}/data/${fmt}-${DATA_VERSION}.tgz)
    set(script
        ${CMAKE_CURRENT_BINARY_DIR}/download_${DATA_PREFIX}_${fmt}.cmake)
    # the marker of the extraction, on which what is made from the instances
    # depends too
    set(DATA_MARKER_${fmt}
        ${CMAKE_CURRENT_SOURCE_DIR}/data/.extracted-${fmt}-${DATA_VERSION})

    file(WRITE ${script} "
file(DOWNLOAD \"${url}\" \"${archive}\" SHOW_PROGRESS STATUS _status)
list(GET _status 0 _code)
list(GET _status 1 _message)
if(NOT _code EQUAL 0)
  # what is on disk is whatever the server answered with, and leaving it there
  # would let the build take it for the archive: the next attempt has to
  # download again rather than fail while extracting it
  file(REMOVE \"${archive}\")
  message(FATAL_ERROR \"cannot download ${url}: \${_message}\")
endif()
file(SIZE \"${archive}\" _size)
if(_size LESS 1024)
  file(REMOVE \"${archive}\")
  message(FATAL_ERROR \"${archive} came down with \${_size} bytes, so the \
server answered but not with the archive: try again later\")
endif()
")

    add_custom_command(
            OUTPUT ${archive}
            COMMAND ${CMAKE_COMMAND} -P ${script}
            COMMENT "Downloading ${fmt}.tgz from the GitLab Package Registry"
            VERBATIM)
    add_custom_target(download_${DATA_PREFIX}_${fmt} DEPENDS ${archive})

    add_custom_command(
            OUTPUT ${DATA_MARKER_${fmt}}
            COMMAND ${CMAKE_COMMAND} -E tar xzf ${archive}
            COMMAND ${CMAKE_COMMAND} -E touch ${DATA_MARKER_${fmt}}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
            COMMENT "Extracting ${fmt}.tgz"
            VERBATIM
            DEPENDS ${archive})
    add_custom_target(extract_${DATA_PREFIX}_${fmt}
                      DEPENDS ${DATA_MARKER_${fmt}})
    # the archive is the output of a command in both targets: it has to be
    # downloaded once, before it is extracted
    add_dependencies(extract_${DATA_PREFIX}_${fmt}
                     download_${DATA_PREFIX}_${fmt})
endforeach ()

# A handful of the textual instances written as netCDF files, which is what a
# tester that reads a Block from a file needs: tools/batch says which ones and
# why, the tests of the module building their own instances instead.
if (UNIX)
    add_custom_command(
            OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated
            COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tools/batch
                    ${CMAKE_CURRENT_BINARY_DIR}/tools/bk2nc4
            COMMAND ${CMAKE_COMMAND} -E touch
                    ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
            COMMENT "Running bk2nc4"
            VERBATIM
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tools/bk2nc4
                    ${DATA_MARKER_txt})
elseif (WIN32)
    add_custom_command(
            OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated
            # using Git for Windows that includes Git Bash
            COMMAND "C:/Program Files/Git/bin/bash.exe"
            "${CMAKE_CURRENT_SOURCE_DIR}/tools/batch" "$<TARGET_FILE:bk2nc4>"
            COMMAND ${CMAKE_COMMAND} -E touch
                    ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
            COMMENT "Running bk2nc4"
            VERBATIM
            DEPENDS $<TARGET_FILE:bk2nc4> ${DATA_MARKER_txt})
endif ()

add_custom_target(run_bk2nc4
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data/nc4/.generated)

add_dependencies(run_bk2nc4 bk2nc4 extract_bk_txt)

# Hard-instance benchmark sources, curated on demand into data/ (gitignored).
# data/fetch-jooken pulls a curated subset of the Jooken-Leyman-De Causmaecker
# hard 0-1 knapsack instances (the .csv carry the published optimum z, read by
# the test in -C mode exactly like the Pisinger ones); data/fetch-txt-large
# the large uncorrelated / subset-sum Pisinger classes. These are the only
# dependency of the test fixtures (see tests/BinaryKnapsackBlock/batch-jooken).
add_custom_target(fetch_jooken_data
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/data/fetch-jooken
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
        COMMENT "Fetching curated Jooken hard knapsack instances")

add_custom_target(fetch_txt_large_data
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/data/fetch-txt-large
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
        COMMENT "Fetching large Pisinger knapsack instances")

# --------------------------------------------------------------------------- #
# Remove from the search paths the cmake directory we added at the beginning.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# --------------------------------------------------------------------------- #
