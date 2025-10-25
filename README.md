# Challenge_2_NLA
2nd NLA Challenge Polimi

Actually with the last request they ask us to comment results so I guess we should write a complete README.md file in order to complete the challenge.

## Point 1

Ag is the adjacency matrix of the graph shown in the Challenge paper. The function `build_adjacency_matrix` builds a Sparse Matrix from a vector of Tuples passed as an argument. 

The Frobenius norm of the matrix $A_g$ is: `4.69042`

## Point 2

After building the $D_g$ and $L_g$ matrix as requested we perform the matrix vector multiplication such that

$$
    y = L_g \cdot x
$$

The Euclidean norm of y: `0`

This result is obvious as we are actually taking the sum of all the terms that are on the same row in Lg, which is the Laplacian.

Lg is SsPD as we get it from the sum of two symmetrix matrices ($D_g - A_g$). (Note that the matrix $L_g$ is semi-positive definite as it presents a 0 eigenvalue. See [point 4](#point-4) for more information.)

## Point 3

Third point asks to compute the biggest and smallest eigenvalue of $L_g$ using the `SelfAdjointEigenSolver` from Eigen.

Minimum eigenvalue of $L_g$ is: `0`

Maximum eigenvalue of $L_g$ is: `5.09259`

## Point 4

As it is said in the Challenge paper $L_g$ is a graph Laplacian which has the smallest eigenvalue equal to zero.

The second smallest eigenvalue is called the Fiedler value. Its corresponding eigenvector carries informations about how the graph is clustered. In particular the Fiedler vector associated to the Fiedler value is the following:

```
Fiedler vector vector: 
-0.440128
-0.385007
-0.233451
-0.385007
 0.128136
 0.273771
 0.350833
  0.37886
 0.311992
```

We can identify two clusters by looking at the sign of the entries of the previous vector.

```
The clusters are: 
Cluster 1 (positive components): 5 6 7 8 9 
Cluster 2 (negative components): 1 2 3 4 
```

This is the result expected because looking at the graph's scheme in the Challenge's paper we clearly recognize two blocks of nodes.

## Point 5

Pretty straight forward point here, you just load the matrix as we always did and report the Frobenius norm here.

Frobenius norm of the matrix: `93.819`

## Point 6

Just repeat the procedure we used for point 2 to obtain the graph Laplacian of the matrix As. In this case the reported output is the following:

$L_s$ is SPD?: No

$L_s$ is symmetric? : Yes

$L_s$ Nonzeros entries = 9153

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

Computed biggest eigenvalue of *Ls.mtx*: `6.013370e+01`

Number of iterations for the method: `2007`

## Point 8
<<<<<<< HEAD
A shift can be added when invoking the previous solver in order to speed up the convergence: 

I made the tests taking 6.01 as a shift, which makes sense because it shifts the maximum eigenvalue closer to zero, and applying the inverse power method is less expensive in terms of iterations. 
This is the command I used:
```
mpirun -n 4 ./eigen1 Ls.mtx eigvec.txt hist.txt -e ii -etol 1e-08 -emaxiter 10000 -shift 6.01
```
And the result was this:
=======
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

Ok it could be me but I guess using the shift meant that you applied the inverse power method. I made the same tests taking 6.01 as a shift, which makes sense because it shifts the maximum eigenvalue closer to zero, and applying the inverse power method is less expensive.
Btw this is the command I used:

```
mpirun -n 4 ./eigen1 Ls.mtx eigvec.txt hist.txt -e ii -etol 1e-08 -emaxiter 10000 -shift 6.01
```

And the result was like this:

>>>>>>> 76ac59ab463d719c50b20627b67fd02eb4554d0f
```
number of processes = 4
matrix size = 351 x 351 (9153 nonzero entries)

initial vector x      : all components set to 1
precision             : double
eigensolver           : Inverse
convergence condition : ||lx-(B^-1)Ax||_2 <= 1.0e-08 * ||lx||_2
matrix storage format : CSR
shift                 : 6.010000e+00
linear solver         : BiCG
preconditioner        : none
eigensolver status    : normal end

Inverse: mode number          = 0
Inverse: eigenvalue           = 6.011940e+00
Inverse: number of iterations = 7
Inverse: elapsed time         = 2.579194e-02 sec.
Inverse:   preconditioner     = 2.639620e-04 sec.
Inverse:     matrix creation  = 2.000000e-07 sec.
Inverse:   linear solver      = 2.540852e-02 sec.
Inverse: relative residual    = 5.457944e-10
```
<<<<<<< HEAD
Which is still worse because of course we are using the inverse power method (it requires to solve a linear system). Number of iterations is good but that's because of the complete different approach of the method we used. I tried other LIS scripts and this one seems the best.
=======

Which is still worse because of course we are using the inverse power method (it requires to solve a linear system). Number of iterations is good but that's because of the method we used. I tried other LIS scripts and this one seems the best. I guess if we are giving a mu as answer this seems fair.
>>>>>>> 76ac59ab463d719c50b20627b67fd02eb4554d0f

## Point 9

Compile the file *etest5.c* by typing:

```
mpicc -DUSE_MPI -I${mkLisInc} -L${mkLisLib} -llis etest5.c -o eigen2
```

The following command runs the LIS script with the **Shift-and-Invert** method which returns the first eigenvalues near the shift (by default is zero, close to the smallest eigenvalue in our case).

```
mpirun ./eigen2 Ls.mtx evals.mtx eigvecs.mtx res.txt iters.txt -ss 2 -e si -etol 1.0e-10
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
Subspace: elapsed time         = 1.996836e-03 sec.
Subspace: number of iterations = 3
Subspace: relative residual    = 4.512717e-12

Subspace: mode number          = 1
Subspace: eigenvalue           = 1.789070e+00
Subspace: elapsed time         = 2.661525e-02 sec.
Subspace: number of iterations = 113
Subspace: relative residual    = 8.965722e-11

eigensolver status    : normal end

Subspace: mode number          = 0
Subspace: eigenvalue           = 5.669404e-04
Subspace: number of iterations = 3
Subspace: elapsed time         = 2.868936e-02 sec.
Subspace:   preconditioner     = 2.136300e-05 sec.
Subspace:     matrix creation  = 1.170000e-07 sec.
Subspace:   linear solver      = 1.922718e-03 sec.
Subspace: relative residual    = 4.512717e-12

```
The second smallest eigenvalue determined is **1.789070e+00** with **113 iterations**

## Point 10

As we already know the dimension of the matrix Ls is 351x351. 

The following command prints the size of a eigenvector corresponding to the second smallest eigenvalue;
```
cat eigvecs.mtx | grep ' 2 ' | wc -l
```

It returns 351 (as expected).

The second command prints the number of negative entries of the previous eigenvector

```
cat eigvecs.mtx | grep ' 2  -' | wc -l
```

`Number of negative entries = 299`.

In order to find the number of positive entries we just need a subtraction

`Number of positive entries = 351 - 299 = 52`.

## Point 11

After constructing the Permutation matrix from the vector of indices relative to the Fiedler's vector, we first defined $A_{ord}$ as:

$$
A_{ord} = P \cdot A_s \cdot P^{T}
$$
Than we counted the numer of non-zero entries in the block of dimension np x np of both As and Aord. The values are reported:

```
Number of nonzero entries in the non-diagonal blocks Aord = 185
Number of nonzero entries in the non-diagonal blocks As = 181
```





