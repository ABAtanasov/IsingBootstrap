# Studying the space of 3D Ising-like CFTs using the conformal bootstrap

This undergraduate research project studies 3D conformal field theories (CFTs) that live in the same "space of theories" as the 3D Ising model's CFT. Specifically, we study theories with two relevant scalars, one being Z_2 even (\epsilon) and the other being Z_2 odd (\sigma). We employ the crossing symmetry constraints on various correlators associated with a given candidate CFT to rule out the space of possible scaling dimensions for these two fields. 

We expand on the earlier work of D. Poland, D. Simmons-Duffin et al., making use of mixed correlator constraints, as well as a "theta scan" method to examine the possible 3-pt coefficient ratios of a given theory. 

This project makes use of T. Ohstuki's cboot module for sage to build the relevant conformal blocks, and D. Simmons-Duffin's sdpb to employ the semidefinite programming for determining the feasibility of a given set of correlators satisfying the boostrap constraints. 
