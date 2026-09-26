# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `BinaryKnapsackSolution::is_dual_feasible()` returns false, the Solution
  holding no dual values [see `Solution::is_dual_feasible()`]

- `bk2nc4`, which writes as a netCDF file a textual instance in the
  Pisinger/Jooken benchmark format or in the native one, so that whoever
  reads a `Block` rather than a knapsack has one to read; it can reverse the
  sense while doing so, changing the sign of the profits, the problem being
  the same one, which is what a component of a decomposition of a
  minimization problem has to be. `tools/batch` writes a handful of the
  curated instances this way and the `run_bk2nc4` target drives it

- `is_sol_feasible()` reads the solution out of the `BinaryKnapsackSolution`
  and checks it against the data of the problem, i.e. the bounds and the
  integrality of the items, the ones that are fixed and the capacity, so that
  the Variable of the `BinaryKnapsackBlock` are neither needed nor touched,
  which `is_sol_feasible_physical()` says; the check that the base class does,
  writing the solution in the Variable and putting back what was there, is
  what a Function revalidating a global pool of hundreds of entries used to
  pay for each of them

### Changed

- the data archive is downloaded by version: `DATA_VERSION` in CMakeLists.txt
  names the version of the Package Registry to read, and the archive and the
  marker of its extraction carry it in their name, so that a tree holding
  an older extraction (the cache of the CI, or a clone extracted before)
  downloads and extracts again instead of running on the old data;
  data/upload-txt publishes the archive under that version

- the items fixed by one call travel as one `GroupModification`, one per call
  and not one loose `Modification` per item, so that a Solver able to write a
  whole set of fixings in one operation can do so, while one that is not
  takes the group apart and sees exactly what it saw before

- the makefile asks for `-O3 -DNDEBUG` and nothing else, the macro of the
  patch for `boost::any` on macOS having no reason to be there since there is
  no `boost::any` left in the core

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

- `chg_weights()`, `chg_profits()` and `chg_capacity()` take their data as a
  `std::span< const double >`, whose length they check against the Range or
  the Subset instead of reading past the end, and are registered in the
  methods factory in that form too; the forms taking an iterator stay, and
  defer to the span ones

- `fix_x()` and `unfix_x()` issue their "abstract" Modification inside a
  GroupModification, one per call, rather than one loose Modification per
  item: a Solver able to execute a whole set of fixings in one operation can
  then do so, while one that is not takes the group apart and sees exactly
  what it saw before

### Fixed

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

- the step that fetches the data archive of this module says what went wrong
  when it goes wrong: the download is checked, an archive that did not arrive
  is removed instead of being left on disk for the build to take for the real
  one, and the message names the URL. A server that answers with an error page
  used to leave a file of a few bytes there, which made the next build fail
  while extracting it, with the message of `tar` and no mention of the
  download

- the sense of the objective survives a netCDF round trip: `serialize()`
  writes the attribute `Sense` when the problem is a minimization one, and
  `deserialize()` reads it, a file without it describing a maximization
  problem, which is what every file written before this did. The sense used
  to be left at its default by `deserialize()`, so that a minimization
  knapsack read back from a file was a maximization one

- `serialize()` wrote the integrality of the items out of a vector it had
  never sized, i.e. past its end, whenever some item was continuous

## [0.4.0] - 2026-09-12

### Added

- `BinaryKnapsackSolver::get_Solution()` and
  `DPBinaryKnapsackSolver::get_Solution()`, which build the
  BinaryKnapsackSolution out of the data structures of the Solver, without
  writing anything into the BinaryKnapsackBlock and therefore without
  requiring any Variable to exist

- accessor and setter to the variables of `BinaryKnapsackSolution`

- `CoreDPBinaryKnapsackSolver`, the dynamic programming Solver enumerating the
  core of the items around the break item, with dominance among the states

- `bk2nc4` and `knap2nc4`, the converters of the Pisinger and of the knapsack
  text instances to netCDF

### Changed

- `GreedyRelaxationBinaryKnapsackSolver` and
  `IncrementalGreedyRelaxationBinaryKnapsackSolver` are ChangeSolver traits,
  with no Block of their own, as the Branch-and-X of the core wants them;
  GroupChange comes from Change.h

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Fixed

- the efficiency order of the items was built with a comparator that a
  negative weight turns into a non-ordering; it is built in one place only

- the package configuration file finds OpenMP, which the library links

## [0.3.2] - 2024-02-27

### Changed

- makefiles and Cmake files updated to new global SMS++ scheme

- documentation updated accordingly

## [0.3.1] - 2022-08-25

### Added

- set_dual() method to BinaryKnapsackBlock

## [0.3.0] - 2022-06-28

### Added

- added BinaryKnapsackBlockChange

- added Physical Representation for fixed variables

### Changed

- default value for reopt parameter of DPBinaryKnapsackSolver

- avoided un-necessary access to BinaryKnapsackBlock by DPBinaryKnapsackSolver

- significant redesign of DPBinaryKnapsackBlock::compute()

### Fixed

- some bugs

## [0.2.0]  - 2021-07-12

### Added

- First complete release with support for the mixed-integer case

### Fixed

- Correctly dealing with negative weights and profits

- Many tests and fixes

## [0.1.0]  - 2021-04-29

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.4.0...develop
[0.4.0]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.3.2...0.4.0
[0.3.2]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.3.1...0.3.2
[0.3.1]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.3.0...0.3.1
[0.3.0]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/binaryknapsackblock/-/tags/0.1.0
