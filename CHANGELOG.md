# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.3.1...develop
[0.3.1]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.3.0...0.3.1
[0.3.0]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/binaryknapsackblock/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/binaryknapsackblock/-/tags/0.1.0
