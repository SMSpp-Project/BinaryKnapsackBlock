# BinaryKnapsackBlock

This project covers two conceptually different things (which may one day be
split to two different projects):

- `BinaryKnapsackBlock`, a SMS++ :Block for the (Mixed) Binary Knapsack
  Problem

- a family of SMS++ :Solver for BinaryKnapsackBlock:

  - `DPBinaryKnapsackSolver`, the standard full-table Dynamic Programming
    approach, with reoptimization support;

  - `ParallelDPBinaryKnapsackSolver`, a DPBinaryKnapsackSolver whose inner
    (capacity) loop can run on parallel engines (OpenMP / FastFlow /
    std::thread); opt-in and serial by default, since this fine-grained
    parallelism does not pay off in practice;

  - `CoreDPBinaryKnapsackSolver`, an exact core / non-dominated-state solver
    in the MINKNAP / COMBO / RECORD line (break-item enumeration with
    dominance, LP / Martello-Toth / surrogate bounds, divisibility and
    aggregation reductions, primal heuristics), orders of magnitude faster
    than the full-table DP;

  - `GreedyRelaxationBinaryKnapsackSolver`, the exact greedy (Dantzig)
    solver of the continuous relaxation, providing true lower / upper bounds
    on the integer optimum and branching on the critical item;

  the last two share the abstract base class `BinaryKnapsackSolver` (raw
  instance mirror with incremental Modification processing, normalized core,
  continuous relaxation).


## Getting started

These instructions will let you build the `BinaryKnapsackBlock` module on
your system.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(BinaryKnapsackBlock)
target_link_libraries(<my_target> SMS++::BinaryKnapsackBlock)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `BinaryKnapsackBlock` and `DPBinaryKnapsackSolver`,
such as the [tester for `BinaryKnapsackBlock` and `DPBinaryKnapsackSolver`
comparing it with
`MILPSolver`](https://gitlab.com/smspp/tests/-/blob/develop/BinaryKnapsackBlock/test.cpp?ref_type=heads),
has to include a "main makefile" of the module, which typically is either
[makefile-c](makefile-c) including all necessary libraries comprised the
"core SMS++" one, or [makefile-s](makefile-s) including all necessary
libraries but not the "core SMS++" one (for the common case in which this is
used together with other modules that already include them). These in turn
recursively include all the required other makefiles, hence one should only
need to edit the "main makefile" for compilation type (C++ compiler and its
options) and it all should be good to go. In case some of the external
libraries are not at their default location, it should only be necessary to
create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Tools

Under `tools/` the module ships `bkbench`, an in-process benchmark driver that
runs any of the solvers over the large public Pisinger / Jooken benchmark
archives (see the Data section), cross-checking each optimum against the
certified one recorded in the data. It defaults to the efficient
`CoreDPBinaryKnapsackSolver`; a comma-separated list selects / cross-checks
other solvers.

You can run it from the `<build-dir>/tools` directory, install it with the
library (see above), or just go in the `tools/` folder and run `make` there
(provided makefiles are properly set, see above). Run it without arguments for
info on its usage:

```sh
bkbench
```


## Data

A small sample of knapsack instances ships with the repo, in the Pisinger text
format ([data/txt](data/txt)); the raw `.csv` kept alongside them embed, per
instance, the certified optimum used by the test suite for the
solver-vs-reference cross-check (see
`tests/BinaryKnapsackBlock/batches/batch-pisinger`).

The full public benchmark archives used by `bkbench` are large and are not
included; populate the gitignored [data/benchmarks](data/benchmarks) folder
via

```sh
cd data
./fetch-benchmarks pisinger jooken
```

which downloads David Pisinger's instances
(https://hjemmesider.diku.dk/~pisinger/codes.html) and the
Jooken-Leyman-Causmaecker instances
(https://github.com/JorikJooken/knapsackProblemInstances).


## Getting help

If you need support, you want to submit bugs or propose a new feature, you
can [open a new
issue](https://gitlab.com/smspp/binaryknapsackblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

## Authors

### Current Lead Authors

- **Federica Di Pasquale**  
  Dipartimento di Informatica  
  Università di Pisa

- **Francesco Demelas**  
  Laboratoire d'Informatique de Paris Nord  
  Universite' Sorbonne Paris Nord

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
