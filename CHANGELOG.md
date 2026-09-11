# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

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

- First complete release with support for the mixed-integer case

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
