# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- a tester of the module in `test/`, run by `ctest -L BinaryKnapsackBlock`,
  which needs nothing but the core: on small instances built in memory, whose
  optimum it finds by enumeration, every Solver of the module has to find the
  optimum (or the value of the relaxation and bounds on the two sides of the
  optimum) with a feasible Solution worth it, over the corners of the data,
  the fixings, the two senses, a sequence of changes made both on the Block
  and through its abstract representation, the netCDF round trip, the R3
  copy and the branching of the relaxation Solvers; the CI of the module
  builds it alone with the core, as a user of the module would

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

- `intReopt` 4 of `CoreDPBinaryKnapsackSolver`: as 3, and moreover an item
  that the Lagrangian bound does not clear is tried again with the
  Martello-Toth bound with that item flipped (the larger of the two
  continuous bounds with the critical item fixed to 0 and to 1), which sees
  that an item as large as the capacity is taken whole or not at all, as in
  the subproblems of the Lagrangian relaxation of facility location
  problems; `get_reopt_outcome()` returns 4 when it was needed

- the certificate of `CoreDPBinaryKnapsackSolver` (`intReopt` 2 and 3)
  covers the changes of the weights too: the items of the previous solution
  may gain weight while it still fits, and the other items may gain weight
  or lose it, an item losing weight being then checked by the Lagrangian
  bound as one whose profit moved

- `CoreDPBinaryKnapsackSolver` reoptimizes a sequence of solves, as chosen by
  `intReopt`: with 1 the previous optimal solution, repaired to the new data,
  is the incumbent the enumeration starts from; with 2 it is moreover
  returned with no solve at all when the changes since cannot have made it
  suboptimal (same core items, a capacity no larger that it still fits, the
  taken ones with the same weight and no less profit, the others with no
  less weight and no more profit); with 3 also when some profit moved against it, provided that for
  each such item a Lagrangian bound with that item flipped cannot beat it,
  the test being skipped for a doubling number of solves after repeated
  failures

- `CoreDPBinaryKnapsackSolver::get_reopt_outcome()` says how much of the
  previous solve the last one reused: nothing, the previous solution as the
  incumbent, or the previous solution itself, returned with no solve because
  the changes could not affect it or because a Lagrangian bound showed it

- the components of the enumeration of `CoreDPBinaryKnapsackSolver` are
  switched by int parameters instead of macros, so that the choice can be
  made per instance in a configuration file: `intSurrogate`,
  `intSurrTrigger` and `intSurrAdapt` for the surrogate relaxation with a
  cardinality constraint and its adaptive rules, `intDPExtension`,
  `intDominanceFix`, `intReductionFix`, `intPrimalHeur` and `intLazyCore`;
  all of them are exact, changing the running time and never the optimum

- `RECORDBinaryKnapsackSolver` and `COMBOBinaryKnapsackSolver` hand the
  integer core of the instance to RECORD and to COMBO, everything else
  (the mirror of the data, its reduction, the reoptimization of `intReopt`
  and the case with continuous variables) being that of
  `CoreDPBinaryKnapsackSolver`, from which they derive. They are only built,
  by CMake and by the makefiles, if the sources are found: `RECORD_ROOT` is a
  checkout of RECORD, which `INSTALL.sh` makes, and `COMBO_ROOT` the
  directory with `combo.c` and `combo.h` of its authors, which are for
  academic or non-commercial use only and are not distributed; each is
  looked for in the variable, in the environment and in the default path of
  `extlib`. The sources are used as they are, but for the `main()` of
  RECORD, which is cut away, and the products of profits and weights in
  COMBO, which are widened to 128-bit integers so that its tests stay exact
  past 2^53; the tester of the module runs them when the build has them

### Changed

- the enumeration of `CoreDPBinaryKnapsackSolver` does less work per state
  and generates fewer states: the break item is found by quickselect, the
  states carry a 64-bit mask of the last items with checkpoints from which
  the solution is rebuilt, the profits are integers whenever the data allow
  it, the core is first solved on a small window around the break item as a
  heuristic, and the items are reduced against the incumbent and fixed by
  dominance while the enumeration grows

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

- `CoreDPBinaryKnapsackSolver` could return the optimum less one: a bound
  computed in floating point that is mathematically an integer could land
  just below it and lose a unit when rounded, which bounds now do past a
  relative tolerance in the direction that can only loosen them; the case is
  in the tester

- `CoreDPBinaryKnapsackSolver::compute()` releases the lock of the Solver
  before passing on an exception thrown by the solve of the core (as the
  external solvers of the derived classes may do), instead of leaving it
  locked

- the data archive is extracted by `cmake -E tar`, which also works with the
  tar of macOS, where the option `--warning=no-unknown-keyword` of GNU tar
  stopped the build.
- the items fixed in `load()` are fixed in their ColVariable too, at the
  value they are fixed to: the ColVariable used to be left free and at 0, so
  that `DPBinaryKnapsackSolver`, which reads the fixings there, took an item
  fixed to 1 for one fixed to 0, and `is_empty()` did not see the fixings at
  all; an R3 copy, which is made by `load()`, had the same defect

- the conditional bounds of the Block, i.e. `get_valid_lower_bound()` and
  `get_valid_upper_bound()`, take each fixed item at the value it is fixed
  to, instead of as a free one: a profitable item of negative weight fixed to
  0 used to be counted in the lower bound, which could then be above the
  optimum

- `fix_x()` on a Subset that is not ordered gives each item the value in the
  same position, as documented: the Subset used to be sorted without the
  values, which then went to other items

- `ParallelDPBinaryKnapsackSolver.h` can be included without the headers of
  FastFlow, the constructor being defined in the .cpp as the destructor is,
  since it destroys the `ff::ParallelFor` if anything throws

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
