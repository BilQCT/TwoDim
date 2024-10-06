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

# TwoDim Background
Simplicial distributions model distributions on a space of outcomes parametrized by a space of measurements. In this code we restrict to two-dimensional measurement spaces and $\mathbb{Z}_d=\{0,1,\cdots,d-1\}$ as the outcomes. In this theory measurement and outcome spaces are modeled by simplicial sets. For an introduction to two-dimensional simplicial sets see [[Okay 2023]](https://arxiv.org/abs/2312.15794).

 
A simplicial set $X$ consists of the following data:

* a collection of sets $X_0,\cdots, X_n, \cdots$
* face maps $d_i:X_n\to X_{n-1}$
* degeneracy maps $s_j:X_{n}\to X_{n+1}$ 

where the face and the degeneracy maps satisfy the simplicial relations. The set $X_n$ represents the set of $n$-simplices. The face and degeneracy maps encode the gluing and collapsing data. A simplex is called non-degenerate if it is not in the image of a degeneracy map.

A simplicial map $f:X\to Y$ consists of functions $f_n:X_n\to Y_n$, $n\geq 0$, compatible with the simplicial structure. Writing $f_\sigma = f_n(\sigma)$ this means that

$ d_i f_\sigma = f_{d_i\sigma} \;\;\;\; s_j f_\sigma = f_{s_j \sigma} $

for all face and degeneracy maps. A simplicial map is determined by $f_\sigma$'s where $\sigma$ is a non-degenerate simplex of $X$.
 

A two-dimensional simplicial set comes with

*  the sets $X_0,X_1,X_2$ of simplices
* the face maps
$$
d_i: X_1 \to X_0 \;\;\;\; i=0,1
$$
and
$$
d_i:X_2\to X_1 \;\;\;\; i=0,1,2
$$
* the degeneracy map
$$
s_0:X_0\to X_1
$$
and
$$
s_0:X_1\to X_2 \;\;\;\; s_1:X_1\to X_2
$$
Two-dimensional simplicial sets will represent our measurement spaces. 

Outcome spaces will be represented by the nerve space $Y=N\mathbb{Z}_d$. We are only concerned with the two-dimensional part of this space:

*  the sets 
$$Y_0=\{\ast\} \;\;\;\; Y_1=\mathbb{Z}_d \;\;\;\; Y_2=\mathbb{Z}_d^2 $$
of simplices
* the face maps
$$
d_0(a,b) = b\;\;\;\; d_1(a,b)=a+b\;\;\;\;d_2(a,b)=a
$$ 
* the degeneracy map
$$
s_0(\ast) = 0.
$$

 

A two-dimensional simplicial distribution on $(X,Y)$ is a simplicial map $p:X\to DY$ where $DY$ is the space of distributions constructed from the outcome space. An $n$-simplex of this space is a distribution on the set of $n$-simplices of the outcome space. For our case a simplicial distribution consists of

* distributions 
$$
p_\sigma = \{p_\sigma^{ab}\}_{a,b\in \mathbb{Z}_d} 
$$
on the set of pairs $(a,b)$ parametrized by the two-simplices $\sigma\in X_2$ satisfying the simplicial relations coming from the face and the degeneracy maps. 

For example, let our measurement space be the diamond space $D$ obtained by gluing two triangles ($2$-simplices) $\sigma$ and $\sigma'$ along a common face, i.e., say $d_1\sigma=d_1\sigma'$. Then a simplicial distribution on $D$ consists of two distributions $p_\sigma$ and $p_{\sigma'}$ such that the marginals along the $d_1$ faces match. 


## Convex polytope of simplicial distributions


The goal is to construct the polytope $P_X$ of simplicial distributions defined on $(X,Y)$ where $X$ is two-dimensional and $Y$ is the nerve space. 
<!--for a simplicial scenario $(X,N\mathbb{Z}_2)$ whose maximal non-degenerate simplices are $2$-dimensional.--> For simplicity, let us consider a two-dimensional space $X$ in which all non-degenerate $1$-simplices $\tau_i$ are the face of some triangle $\sigma$ so that $\tau_i = d_j\sigma$. This is not particularly restrictive as it is a basic fact about simplicial distributions that simplicial distributions restricted to one-dimensional spaces yield only "classical" distributions, which are not particularly interesting.

By the Minkowski-Weyl theorem a polytope can be defined as the intersection of a finite number of half-space inequalities or as the convex hull of a finite number of extreme points. The former is called the $H$-representation while the latter is called $V$-representation of the polytope. The problem of converting from $H$ to $V$ representation is called the vertex (facet) enumeration problem. The polytope $P_X$ can be defined in its $H$-representation as follows. For $(a,b)\in \mathbb{Z}_2^2$ and non-degenerate $2$-simplex $\sigma$ let us consider the distribution $p_\sigma$ on each triangle. Since we work with the semi-ring $\mathbb{R}_{\geq 0}$ we require that $p_\sigma^{ab}\geq 0$ for all $(a,b)\in \mathbb{Z}_2^2$ and non-degenerate $2$-simplex $\sigma$. This yields an unbounded polyhedron (a cone) given by the positive orthant of a Euclidean space. Moreover, each distribution $p_\sigma^{ab}$ is normalized $\sum_{a,b} p_\sigma^{ab} =1 $ and there are additional compatibility conditions $d_i \sigma = d_j\sigma^\prime$ leading to the constraints

$$\sum_{ab\in\mathbb{Z}_2^2~:~d_i(ab)=c} p_\sigma^{ab} = \sum_{ab\in\mathbb{Z}_2^2~:~d_j(ab)=c} p_{\sigma^\prime}^{ab}.$$


The intersection of the cone with these affine subspaces yields a polytope $P_X$.

### Degeneracies

<!--A crucial distinction between simplicial sets and abstract simplicial complexes is the notion of degeneracies: $s_j:X_n\to X_{n+1}$. If a simplex $\sigma_n \in X_n$ is in the image of a degeneracy map acting on a non-degenerate simplex $\sigma_{n-1}$, i.e. $\sigma = s_j(\sigma_{n-1})$, then this degenerate simplex is not pictured in the geometric realization of $X$. Nonetheless, degeneracies become quite important when we wish to perform natural topological operations, such as collapsing. Consider a collapse map

$$\pi: \Delta^1\to\Delta^0.$$

In the framework of simplicial complexes we could not implement this collapse as a simplicial map because there is no notion of mapping a $1$-simplex to a $0$-simplex. Degeneracies solve this issue since the simplicial set for $\Delta^0$ includes $n$-simplices at every level $(\Delta^0)_n$ that are just images of the non-degenerate point $c\in (\Delta^0)_0$ under the degeneracy map. By mapping a non-degenerate $1$-simplex $\tau\in (\Delta^1)_1$ to a degenerate $1$-simplex $s_j(c)\in (\Delta^0)_1$ we implement the desired collapse since degenerate simplices do not appear in the geometrical realization of $\Delta^0$.-->

We consider the behavior of degeneracies in the case of two-dimensional simplicial distributions. <!--Let $p:X\to DY$ be a simplicial distributoin where $X$ is two-dimensional.--> Suppose we have a non-degenerate $2$-simplex $\sigma\in X_2$ such that one of its faces is a degenerate $1$-simplex, i.e. $d_i(\sigma) = s_0(c)$ for some $c\in X_0$. Since $p$ is a simplicial map we also have

$$d_i(p_\sigma) = s_0(p_c).$$

<!--For an outcome space $N\mathbb{Z}_d$ we can associate $p(\sigma)$ with a tuple $(p^{ab})$ of length $d\times d$. There are three face maps given by

$$ \left(d_i(p_\sigma)\right)(r) =%
\begin{cases}
      \sum_a p_\sigma(ar) & i=0 \\
      \sum_{a+b = r ~\text{mod}~d}p_\sigma(ab) & i =1 \\
      \sum_a p_\sigma(ra) & i = 2~.
    \end{cases}$$-->

Conversely we have that there is a unique distribution on a $0$-simplex given by $p_c = 1$. The degeneracy map acting on this distribution is given by a tuple $s_0(p_c) = \delta^0 $, which is a delta-distribution such that $\delta^0(a) = 1$ for $a = 0$ and $\delta^0(a) =0$ otherwise. Thus the relation $d_i(p_\sigma) = s_0(p_c)$ yields

\begin{align}
\sum_a p_\sigma(a0) &=& 1 \quad\quad (i = 0),\notag\\
\sum_{a+b = 0 ~\text{mod}~d}p_\sigma(ab) &=& 1 \quad\quad (i = 1),\notag\\
\sum_a p_\sigma(0a) &=& 1 \quad\quad (i = 2).\notag\\
\end{align}

Due to normalization of $p_\sigma(ab)$ this implies that all probabilities not appearing in the above expression should be set to zero in the case that a face of $\sigma$ is degenerate. This is imposed by introducing additional linear equalities into the description of the polytope $P_X$.


## Edge coordinates

<!--For simplicial scenarios the map from the boundary of a triangle to the triangle is injective. In particular, for a space $X$ consisting of a single non-degenerate $2$-simplex $\sigma$, consider a simplicial distribution $p:X\to DN\mathbb{Z}_2$, which we denote $p^{ab}_\sigma := p(\sigma)(ab)$.--> The boundary $\partial\sigma$ consists of three edges $x_0,x_1,x_2$ which are the faces $d_0\sigma,d_1\sigma,d_2\sigma$, respectively. Distributions $p_{x_i}^a := p(x_i)(a)$ on the boundary are obtained by marginalization. More explicitly,

\begin{equation}
p_{x_i}^c=\sum_{a,b\in \mathbb{Z}_2~: ~d_i(ab)=c}p^{ab}_\sigma~.\notag
\end{equation}

### Edge coordinates: marginals:

Recall that the marginal distributions $p_{x_i}^a$ are normalized so that $p_{x_i}^0+p_{x_i}^1=1$. <!-- and, without loss of generality,--> The distribution $p_{x_i}$ can be characterized by a single parameter $p_{x_i}^0$ since $p_{x_i}^1=1-p_{d_i\sigma}^0$. Similarly, the distribution $p_\sigma^{ab}$ is itself normalized, thus $p_\sigma^{ab}$ can be characterized by three independent parameters <!--, which we can take to be--> $p_{x_i}^0$ where $i=0,1,2$. To see this, using the formula above, note that

\begin{align}
p_{x_0}^0 &=& p^{00}_\sigma+p^{10}_\sigma,\\
p_{x_1}^0 &=& p^{00}_\sigma+p^{11}_\sigma,\\
p_{x_2}^0 &=& p^{00}_\sigma+p^{01}_\sigma.
\end{align}

Together with normalization we obtain a description of $p_\sigma^{ab}$ in terms of its marginal distributions

$$ p_\sigma^{ab} = \frac{1}{2}\left ( p_{x_0}^a + p_{x_1}^{a+b+1} + p_{x_2}^b \right).$$

### Polytope of simplicial distributions

The distribution $p_\sigma^{ab}$ therefore is uniquely determined by its edge marginals $p_{x_i}^0$, and is automatically normalized. Moreover, <!--for two-dimensional simplicial scenarios $(X,N\mathbb{Z}_2)$,--> non-trivial compatibility (or nonsignaling) constraints arise only when simplices are glued along edges. (Gluing along vertices is possible but does not induce any compatibility constraints.) Thus compatibility of simplicial distributions is guaranteed with edge coordinates. In other words, for two-dimensional spaces with outcomes in $\mathbb{Z}_d$, simplicial distributions are completely characterized by the marginal distributions of the edges.

Let us denote $X_n^\circ$ to be the set of non-degenerate $n$-simplices in $X$. (For our purposes this is only non-empty for $n\leq 2$.) Rather than letting the polytope $P_X$ reside in an ambient space $\mathbb{R}^{|\mathbb{Z}_d^2|\times |X_2^\circ|}$, with coordinates given by the $p_\sigma^{ab}$ $(\sigma\in X_2^\circ)$, it can instead be embedded in a real Euclidean space $\mathbb{R}^{|X_1^\circ|}$ with coordinates $p_\tau^0$ $(\tau\in X_1^\circ)$ where $|X_1^\circ|\leq |\mathbb{Z}_d^2|\times |X_2^\circ|$. (When a polytope $P$ with affine span of dimension $N$ is embedded in $\mathbb{R}^N$, it is called full-dimensional.) The polytope $P_X$ is given in its $H$-representation by the following set of inequalities:

$$ P_X = \{x\in \mathbb{R}^{|X_1^\circ|}~:~p_\sigma^{ab}\geq 0\}$$

where $p_\sigma^{ab}$ is given by the formula above.

### Edge coordinates: expectations:

Some computations are made simpler by a change of coordinates. We interpret the edges $x_i$ as measurements from which we obtain outcomes $(-1)^a$ ($a\in\mathbb{Z}_d$) that occur with probability $p_{x_i}^a$. The expected value (or expectation) is then given by

$$\bar{x}_i := \sum_a (-1)^a p_{x_i}^a = p_{x_i}^0 - p_{x_i}^1.$$

Using the normalization of $p_{x_i}^a$ we have that $\bar{x}_i=2p_{x_i}^0-1$ yielding $p_{x_i}^0=(\bar{x}_i+1)/2$. Plugging into the expression for $p_\sigma^{ab}$ above, we have that:

\begin{equation}
p_\sigma^{ab} = \frac{1}{4}\left( 1+(-1)^a\bar{x}_0+(-1)^{a+b}\bar{x}_1+(-1)^b\bar{x}_2 \right).\notag
\end{equation}


## Twisted simplicial distributions

Motivated by quantum mechanics, it is possible to have so-called twisted distributions on a simplicial scenario. We again restrict our attention to two-dimensional measurement spaces <!--$X$--> and <!--$\mathbb{Z}_2$--> the nerve space as the outcomes. Let $\sigma$ be a non-degenerate triangle with non-degenerate edges $x_i = d_i\sigma$. Denoting the set of non-degenerate edges on the boundary of $\sigma$ by $\partial\sigma$, we consider a function $\beta:\partial\sigma\to\mathbb{Z}_2$. Simplicial maps $r:X\to Y$ preserve the simplicial structure so that if $r(\sigma) = (a,b)$ then $x_i\mapsto d_i(a,b)$. Allowing for twisting amounts to modifying this relation according to
$$r(x_i) = d_i(a,b)+\beta(x_i).$$
Such twisted distributions appear quite naturally in quantum mechanics, where the appearance of a non-trivial $\beta$ function is an indicator of non-classical behavior, often called quantum contextuality. For more on the connection to quantum mechanics, see [[Okay, et al. 2017]](https://arxiv.org/abs/1701.01888)[[Okay, et al. 2023]](https://quantum-journal.org/papers/q-2023-05-22-1009/).



## Conventions and notation for two-dimensional distributions
We input a two-dimensional simplicial set as:
-  a set $X_1^{\circ}$ of non-degenerate edges $\tau$, which we encode individually as an array: $\tau_i = [i,[d_0\tau,d_1\tau]]$,

- we index **degenerate** edges $s_0(c)$ $(c\in X_0)$ by negative numbers which we formally write as $[-i,[c,c]]$, <!--although we do not need to enumerate these in the definition of a simplicial set.-->

- a set $X_2^{\circ}$ of non-degenerate triangles $\sigma$, each of which we encode as an array: $\sigma_i = [i,[d_0\sigma,d_1\sigma,d_2\sigma]]$.

#### Twisting

We may additionally choose to construct a twisted scenario, in which an edge $\tau$ in a triangle $\sigma$ may have a "twisted" outcome assignment $\tau\mapsto r(\tau)+\beta(\tau)$, as described above. Note that it suffices to twist just one edge. We therefore encode the twisting by an array $\mathcal{T} = [T_i]$ of elements $T_i = [i,\beta,k]$, where $i$ is the identifier for $\sigma_i\in X_2^{\circ}$, $\beta\in \mathbb{Z}_d$ identifies the twisting on the $d_k(\sigma_i)$ edge.

**NOTE:** We use Julia indexing so that for $d_k:X_2\to X_1$ we use $k=1,2,3$ rather than $k=0,1,2$ in the literature.

**NOTE:** By default we assume no twisting.
