**NOTICE: Code usupported with Julia > v0.7**
---

A Julia based code to calculate atomic bond-angle distribution from
molecular dynamics simulations. 

>Copyright (2015) Stefan Bringuier
>This software is distributed under the GNU General Public License.
>See LICENSE file for more details.

Program Hierarchy:
README		this file
USAGE		how to use this code
LICENSE		the GNU General Public License (GPL)
```shell
src/		source files
		--neighbor.jl		calculate neighbor list
		--angle.jl 		calculate bond angles
		--lattice.jl		generate various lattices
		--test.jl 		contains unittest for functions
		--readlammps.jl		read LAMMPS dump and data format 
		--writelammps.jl 	write LAMMPS dump and data format
		--orderparam.jl         measure the order parameter of NiTi from Mutter et. al APL 2011					   
```

This code is a collection of Julia modules which contain various functions related to analyzing molecular dynamics simulation trajectories from the LAMMPS package. The code was written for two purposes:

1. For me to learn Julia which I believe is a very very promising language. 
2. Provide an optimized version with regards to the my python based bond angle analysis code.

Since Julia makes use of JIT compiling I'm hoping that this will provide the required performance increase (without resorting to parallelization) that I desire, especially for the bond angle routine which scales poorly. If anyone has a better algorithm for this please let me know, it would be much appreciated. Additionally, if you think there are better ways to code in Julia in contrast to how I have done them, please let me know.
