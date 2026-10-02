# test

A tester for `BinaryKnapsackBlock` and its Solvers that needs nothing but the
core SMS++ library.

It builds in memory small instances of the (mixed) binary knapsack, whose
optimum it finds by enumerating the 0-1 items and filling the continuous ones
greedily by efficiency for each assignment of the others; the same
enumeration with every item continuous gives the continuous relaxation. The
exact Solvers (`DPBinaryKnapsackSolver`, also with the reoptimization on,
`ParallelDPBinaryKnapsackSolver` with each of its engines and
`CoreDPBinaryKnapsackSolver`) have to find the optimum, with a feasible
`Solution` worth it, or say that the instance is empty; the relaxation ones
(`GreedyRelaxationBinaryKnapsackSolver` and its incremental variant) have to
find the value of the relaxation, bounds on the two sides of the optimum and a
feasible rounded `Solution` worth the bound it gives.

The instances cover random data of both signs and both senses, with
continuous and fixed items, and the corners of the data: capacity 0, no item
fitting, every item fitting, items of weight 0, a negative capacity, an empty
problem and no item at all. On a Block with every Solver attached, a sequence
of changes of profits, weights, capacity, fixings, sense and integrality is
solved again after each change, and so are the same changes made through the
abstract representation (the coefficients of the `LinearFunction` of the
`Objective` and of the constraint, its right-hand side, the fixing of a
`ColVariable`, the sense of the `Objective`). The tester also checks the
conditional bounds and the emptiness of the Block, the abstract
representation against the data, the netCDF round trip through a temporary
file, the R3 copy with the mapping of the solutions, the refusal of weights
that are not integer by the DP Solvers, and the branching of the relaxation
Solvers, whose two children have to hold the optimum and be undone.

The Solvers are attached by the `BlockSolverConfig` files of this directory,
which the tester only reads and applies: `BSPar-relax.txt` the relaxation
ones, `BSPar-exact.txt` the exact ones, `BSPar-RECORD.txt` and
`BSPar-COMBO.txt` the external ones when their sources were given, and
`BSPar-reopt.txt` the core DP with `intReopt` 3, 4 and 5, whose outcomes after
each change are checked against each other; the `ComputeConfig` fragments
`DPCfg.txt`, `PDPCfg.txt` and `ReoptCfg.txt` hold the parameters, which the
variants override. The tester therefore runs from this directory, which is
also the one `ctest` uses.

The exit code is 0 when every check passes, printing `All tests passed!!`, and
1 otherwise. The `makefile` builds the executable including the
`BinaryKnapsackBlock` module and the core SMS++ library.


## Authors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
