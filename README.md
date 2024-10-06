# TwoDim

Simplicial distributions are new mathematical objects introduced in [[Okay, et al. 2023]](https://quantum-journal.org/papers/q-2023-05-22-1009/) that provide a rigorous approach to studying quantum contextuality based on the theory of simplicial sets. In quantum mechanics not all observables can be measured simultaneously, only those that commute. Given a set of measurements, simplicial sets encode these commutivity relations, which are then, in turn, reflected in the measurement statistics using those observables.

This repository contains easy-to-use code for computing geometric and combinatorial properties of the convex set of simplicial distributions when the underlying spaces are two-dimensional in a sense that we define below. We also have the [Examples](Examples.ipynb) Jupyter notebook which demonstrates how to use the code in practice. 

The main working example for two-dimensional distributions, which may generally be twisted, comes from the Mermin square [[Link]](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.65.803). We treat this example in a separate notebook.


# Instructions

To create an environment for running the TwoDim Jupyter notebooks do the following:

In main TwoDim directory enter Julia REPL and use the following commands:

using Pkg
Pkg.generate("MyProject")
Pkg.activate("MyProject")

Pkg.add("YourPackage")

## Requirements
GAP
Polymake
Combinatorics
MacroTools
