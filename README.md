# Challenge_2_NLA
2nd NLA Challenge Polimi

## Request 7
In order to make the graph laplacian matrix invertible, add a small perturbation to the first diagonal entry of $L_s$ , namely $L_s (1, 1) = L_s (1, 1) + 0.2$. Export Ls in the *.mtx* format and move it to the `lis-2.1.10/test` folder. Using the proper iterative solver available in the LIS library compute the largest eigenvalue of Ls up to a tolerance of $10^{−8}$. Report the computed eigenvalue and the iterations counts.

Since the largest eigenvalue is requested, it is possible to apply the power method. This can be achieve by the LIS command: 
```
mpirun -n 4 ./eigen1 Ls_perturbed.mtx eigvec.txt hist.txt -e pi -etol 1.e-8 -emaxiter 3000
```
whose console output is: 
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
Power: elapsed time         = 2.189406e-02 sec.
Power:   preconditioner     = 0.000000e+00 sec.
Power:     matrix creation  = 0.000000e+00 sec.
Power:   linear solver      = 0.000000e+00 sec.
Power: relative residual    = 9.940435e-09

```
