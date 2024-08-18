# The NAS Parallel Benchmarks for evaluating C++ parallel programming frameworks on shared-memory architectures

The NPB's Fortran codes were carefully ported to **C++** and are fully compliant to the **NPB3.4.1** version ([NPB official webpage](https://www.nas.nasa.gov/publications/npb.html)). Our [paper](https://doi.org/10.1016/j.future.2021.07.021) contains abundant information on how the porting was conducted and discusses the outcome performance we obtained with **NPB-CPP** on different machines (Intel Xeon, AMD Epyc, and IBM Power8) and compilers (GCC, ICC, and Clang). Results showed that we achieved similar performance with **NPB-CPP** compared to the original **NPB**. **You can use our paper, along with the official reports, as a guide to assess performance using the NPB suite**.

## How to cite this work
  
[[DOI]](https://doi.org/10.1016/j.future.2021.07.021) J. Löff, D. Griebler, G. Mencagli et al., **The NAS Parallel Benchmarks for evaluating C++ parallel programming frameworks on shared-memory architectures**, *Future Generation Computer Systems (FGCS)* (2021)

*This is a repository aimed at providing parallel codes with different C++ parallel programming APIs for the NAS Parallel Benchmarks (NPB). You can also contribute with this project, writing issues and pull requests.*

The conventions we used in our porting can be found [here](notes-conventions.md)


    ===================================================================
      NAS Parallel Benchmarks in C++ using OpenMP, FastFlow and Intel TBB

        This project was conducted in the Parallel Applications
        Modelling Group (GMAP) at PUCRS - Brazil.

        GMAP Research Group leader:
            Luiz Gustavo Leão Fernandes

        Code contributors: 
            Dalvan Griebler (PUCRS)
            Gabriell Araujo (PUCRS)
            Júnior Löff (PUCRS)

      In case of questions or problems, please send an e-mail to us:	
        dalvan.griebler@acad.pucrs.br
        gabriell.araujo@edu.pucrs.br			
        junior.loff@edu.pucrs.br				

      We would like to thank the following researchers for the 
      fruitful discussions:
          Gabriele Mencagli	(UNIPI)
          Massimo Torquati	(UNIPI)
          Marco Danelutto (UNIPI)
    ===================================================================


### Folders inside the  project:

**NPB-SER** - This directory contains the sequential version.

**NPB-OMP** - This directory contains the parallel version implemented with OpenMP (based in the original NPB version).

**NPB-TBB** - This directory contains the parallel version implemented with Threading Building Blocks.

**NPB-FF** - This directory contains the parallel version implemented with FastFlow.

# The Five Kernels and Three Pseudo-applications

Each directory is independent and contains its own implemented version of the kernels and pseudo-applications:

## Kernels

	EP - Embarrassingly Parallel, floating-point operation capacity
	MG - Multi-Grid, non-local memory accesses, short- and long-distance communication
	CG - Conjugate Gradient, irregular memory accesses and communication
	FT - discrete 3D fast Fourier Transform, intensive long-distance communication
	IS - Integer Sort, integer computation and communication

## Pseudo-applications

	BT - Block Tri-diagonal solver
	SP - Scalar Penta-diagonal solver
	LU - Lower-Upper Gauss-Seidel solver

*Tip: The pseudo-applications' performance is bounded to the sequential partial differential equation (PDE) solver*

# Software Requirements

*Warning: our tests were made with GCC-9 and ICC-19*

# How to Compile 

Enter the directory from the version desired and execute:

`$ make _BENCHMARK CLASS=_WORKLOAD`

_BENCHMARKs are: 
		
	EP, CG, MG, IS, FT, BT, SP and LU 
																										
_WORKLOADs are: 
	
	Class S: small for quick test purposes
	Class W: workstation size (a 90's workstation; now likely too small)	
	Classes A, B, C: standard test problems; ~4X size increase going from one class to the next	
	Classes D, E, F: large test problems; ~16X size increase from each of the previous Classes  


Command example:

`$ make ep CLASS=A`

# How to Execute

Binaries are generated inside the bin folder

Command example:
	
`$ ./bin/ep.A`

# Compiler and Parallel Configurations

Each folder contains a default compiler configuration that can be modified in the `config/make.def` file.
You must use this file if you want to modify the target compiler, flags or links that will be used to compile the applications.

## Parallel Execution

### Using and configuring the used parallel programming frameworks

The repository already has an additional directory `libs` with the FastFlow and Intel TBB libraries.

For TBB you need to compile the library and load the environment variables, therefore, enter `libs/tbb-2020.1` and execute the following command:

`$ make`

This command will generate a folder inside `libs/tbb-2020.1/build`. Finally, you can load TBB vars within the script `tbbvars.sh`, for example, executing the following command in your terminal:

`$ source libs/tbb-2020.1/build/linux_intel64_gcc_cc7.5.0_libc2.27_kernel4.15.0_release/tbbvars.sh`

### Setting the degree of parallelism (NUM_THREADS)

The degree of parallelism can be set using the `*RUNTIME*_NUM_THREADS` environment variable.

Command example:
		
`$ export OMP_NUM_THREADS=32`
or
`TBB_NUM_THREADS` and `FF_NUM_THREADS`
