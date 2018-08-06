# Studying the space of 3D CFTs with a Z_2 Symmetry

*Code authors: Alex Atanasov (alex.atanasov@yale.edu) and Aaron Hillman (aaron.hillman@yale.edu)* 

*Supervised by David Poland (david.poland@yale.edu)*

Package used for obtaining the results in [1807.05702](https://arxiv.org/abs/1807.05702).

This undergraduate research project studies 3D conformal field theories (CFTs) with a parity symmetry. Specifically, we used this package initially for studying theories with two relevant scalars, one being parity even (\epsilon) and the other being parity odd (\sigma), giving theories with the same operator structure as the 3D Ising model. After this, we lowered the gap assumptions on the next lowest-lying parity odd scalar to make the operator relevant, and used this to constrain the scaling dimensions of the supersymmetric 3D Ising model. In both cases, we applied the usual crossing symmetry constraints on the four point functions associated with a given candidate CFT to rule out the space of possible scaling dimensions for the theories.

We expand on the earlier work of D. Poland, D. Simmons-Duffin et al., making use of mixed correlator constraints, as well as a "theta scan" method to examine the possible 3-pt coefficient ratios of a given theory. 

This project makes use of T. Ohstuki's `cboot` module for sage to build the relevant conformal blocks, and D. Simmons-Duffin's `sdpb` to employ the semidefinite programming for determining the feasibility of a given set of correlators satisfying the bootstrap constraints. 

### Installation requirements: 

* The most recent version of [SDPB](https://github.com/davidsd/sdpb) 
installed, with `.sdpb` placed
in the main directory of this project, together with all the
requisite installations for `sdpb` to run successfully
* Version 6.8-7.0 of [SageMath](http://www.sagemath.org/) installed
* [`cboot`](https://github.com/tohtsky/cboot) installed in the
 appropriate Sage module directory