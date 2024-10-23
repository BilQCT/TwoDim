# TwoDim

Simplicial distributions are new mathematical objects introduced in [[Okay, et al. 2023]](https://quantum-journal.org/papers/q-2023-05-22-1009/) that provide a rigorous approach to studying quantum contextuality based on the theory of simplicial sets. In quantum mechanics not all observables can be measured simultaneously, only those that commute. Given a set of measurements, simplicial sets encode these commutitivity relations, which are then, in turn, reflected in the measurement statistics using those observables.

This repository contains easy-to-use code for computing geometric and combinatorial properties of the convex set of simplicial distributions when the underlying spaces are two-dimensional in a sense that is defined in [Examples.ipynb](./Examples.ipynb).  Properties of two-dimensional simplicial distributions have been studied in some detail. In [[Kharoof, et al. 2023]](https://arxiv.org/abs/2306.01459) the Bell inequalities for certain types of two-dimensional scenarios, called flower scenarios, were explicitly characterized, whereas in [[Okay 2023]](https://arxiv.org/abs/2312.15794) graph-theoretical constructions were used to characterize the vertices of two-dimensional twisted scenarios.

The [Examples](Examples.ipynb) Jupyter notebook contains a technical primer on simplicial distributions, with an emphasis on those arising from two-dimensional scenarios, and also demonstrates how to use our code in practice. The main working example for two-dimensional distributions, which may generally be twisted, comes from the Mermin square [[Mermin 1993]](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.65.803). The corresponding polytope of simplicial distributions for the Mermin square scenario was first studied in [[Okay, et al. 2023]](https://arxiv.org/abs/2210.10186). We treat this example in more detail in a separate notebook.


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

**Acknowledgments:** This work is supported by the US Air Force Office of Scientific Research
under award number FA9550-21-1-0002 and the Digital Horizon Europe project FoQaCiA, GA no.101070558.
