# BinaryKnapsackBlock 

This project covers two conceptually different things (which may one day be
split to two different projects):

- `BinaryKnapsackBlock`, a SMS++ :Block for the Binary Knapsack Problem

- `DPBinaryKnapsackSolver`, a SMS++ :Solver for BinaryKnapsackBlock based on
  a trivial implementation of the standard Dynamic Programming approach


## Getting started

These instructions will let you build `BinaryKnapsackBlock` and 
`DPBinaryKnapsackSolver` on your system.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
sudo make install
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


## Getting help

If you need support, you want to submit bugs or propose a new feature, you
can [open a new
issue](https://gitlab.com/smspp/binaryknapsackblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

### Current Lead Authors

- **Federica Di Pasquale**  
  Dipartimento di Informatica  
  Università di Pisa

- **Francesco Demelas**  
  Laboratoire d'Informatique de Paris Nord  
  Universite' Sorbonne Paris Nord

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
