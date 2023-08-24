# BinaryKnapsackBlock 

This project covers two conceptually different things (which may one day be
split to two different projects):

- `BinaryKnapsackBlock`, a SMS++ :Block for the Binary Knapsack Problem

- `DPBinaryKnapsackSolver`, a SMS++ :Solver for BinaryKnapsackBlock based on
  a trivial implementation of the standard Dynamic Programming approach


## Getting started

These instructions will let you build BinaryKnapsackBlock and 
DPBinaryKnapsackSolver on your system.


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
to use CMake. General instructions are:

- The arrangements of folders must be that envisioned by the
  [Umbrella SMS++ Project](https://gitlab.com/smspp/smspp-project)

- The main step is to edit the makefiles into ../extlib/. There is one for
  each of the external libraries that any module requires, starting with

  = [Boost](https://www.boost.org)

  = [Eigen](http://eigen.tuxfamily.org)

  = [netCDF-C++](https://www.unidata.ucar.edu/software/netcdf)

  that are required by the "core" SMS++ library and therefore by everyone.
  Setting the

```make
lib*INC = -I< paths to include files directories >
lib*LIB = -L< paths to lib files directories > -l< libs >
```

  in each allows one to set any non-standard path if the library is not
  installed in the system (or leave them empty if they are).

- The "core" SMS++ classes have a makefile for building the corresponding
  library in

```sh
SMS++/lib/makefile-lib
```

  The makefile allow to choose the compiler name and the optimization/debug.
  This builds the lib/libSMS++.a that can be linked upon. Also, the

```sh
SMS++/lib/makefile-inc
```

  file is provided for allowing external makefiles to ensure that the library
  is up-to-date (useful in case one is actually developing it).

  Note that the "basic" makefile macros

```make
CC =
SW =
```

  for setting the c++ compiler and its options are "automatically forwarded"
  from the makefile to these of the other SMS++ components, and therefore
  (possibly at the cost of a make clean) ensure consistency during the
  building process.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/binaryknapsackblock/-/issues/new).

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
