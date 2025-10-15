# Challenge_2_NLA
2nd NLA Challenge Polimi

Actually with the last request they ask us to comment results so I guess we should write a complete README.md file in order to complete the challenge.

## Point 7

After adding the small perturbation to the value in position (1, 1) and saving the matrix in Matrix Market format, we run the LIS script that computes the biggest eigenvalue of a matrix with a given tolerance. After compiling the source code with the following command.
```
mpicc -DUSE_MPI -I${mkLisInc} -L${mkLisLib} -llis etest1.c -o eigen1
```
We computed the biggest eigenvalue of Ls.mtx using the power method and setting a tolerance of 1e-08:
```
mpirun -n 4 ./eigen1 Ls.mtx eigvec.txt hist.txt -e pi -etol 1e-08 -emaxiter 10000
```
Note that, even if it wasn't requested, the max iteration parameter was set to 10000 because the method didn't reach the stopping condition in the default 1000 iterations set by LIS.
The output of the execution was the following:
```
number of processes = 4
matrix size = 351 x 351 (9153 nonzero entries)

initial vector x      : all components set to 1
precision             : double
eigensolver           : Power
convergence condition : ||lx-(B^-1)Ax||_2 <= 1.0e-08 * ||lx||_2
matrix storage format : CSR
shift                 : 0.000000e+00
eigensolver status    : normal end

Power: mode number          = 0
Power: eigenvalue           = 6.013370e+01
Power: number of iterations = 2007
Power: elapsed time         = 6.756434e-03 sec.
Power:   preconditioner     = 0.000000e+00 sec.
Power:     matrix creation  = 0.000000e+00 sec.
Power:   linear solver      = 0.000000e+00 sec.
Power: relative residual    = 9.940435e-09
```
Computed biggest eigenvalue of Ls.mtx: `6.013370e+01`
Number of iterations for the method: `2007`

## Point 8
A shift can be added when invoking the previous solver in order to speed up the convergence: 
```
mpirun -n 4 ./eigen1 Ls.mtx eigvec.txt hist.txt -e pi -etol 1e-08 -emaxiter 10000 -shift 29.55
```
whose output is: 
```
number of processes = 4
matrix size = 351 x 351 (9153 nonzero entries)

initial vector x      : all components set to 1
precision             : double
eigensolver           : Power
convergence condition : ||lx-(B^-1)Ax||_2 <= 1.0e-08 * ||lx||_2
matrix storage format : CSR
shift                 : 2.955000e+01
eigensolver status    : normal end

Power: mode number          = 0
Power: eigenvalue           = 6.013370e+01
Power: number of iterations = 1063
Power: elapsed time         = 1.176979e-02 sec.
Power:   preconditioner     = 0.000000e+00 sec.
Power:     matrix creation  = 0.000000e+00 sec.
Power:   linear solver      = 0.000000e+00 sec.
Power: relative residual    = 9.972897e-09

```

[NB: I tried various shifts, 29.55 seems to be the best in term of reducing the number of iterations without affecting too much the relative residual (and without causing a change in the computed eigenvalue!).]

## Point 9

This command runs the LIS script and finds the second smallest (positive) eigenvalue of Ls matrix with a 1.e-10 tolerance.

```
mpirun -n 4 ./eigen1 Ls_perturbed_cp.mtx eigvec.txt hist.txt -e si -ss 2 -etol 1e-10
```

with the following result: 

```
number of processes = 4
matrix size = 351 x 351 (9153 nonzero entries)

initial vector x      : all components set to 1
precision             : double
eigensolver           : Subspace
convergence condition : ||lx-(B^-1)Ax||_2 <= 1.0e-10 * ||lx||_2
matrix storage format : CSR
shift                 : 0.000000e+00
inner eigensolver     : Inverse
linear solver         : BiCG
preconditioner        : none
size of subspace      : 2

compute eigenpairs in subspace:

Subspace: mode number          = 0
Subspace: eigenvalue           = 5.669404e-04
Subspace: elapsed time         = 2.074572e-03 sec.
Subspace: number of iterations = 3
Subspace: relative residual    = 4.512717e-12

Subspace: mode number          = 1
Subspace: eigenvalue           = 1.789070e+00
Subspace: elapsed time         = 2.323120e-02 sec.
Subspace: number of iterations = 113
Subspace: relative residual    = 8.965722e-11

eigensolver status    : normal end

Subspace: mode number          = 0
Subspace: eigenvalue           = 5.669404e-04
Subspace: number of iterations = 3
Subspace: elapsed time         = 2.541831e-02 sec.
Subspace:   preconditioner     = 1.465300e-05 sec.
Subspace:     matrix creation  = 8.700000e-08 sec.
Subspace:   linear solver      = 2.006257e-03 sec.
Subspace: relative residual    = 4.512717e-12

```

The second smallest eigenvalue determined is **1.789070e+00** with **113 iterations** 




