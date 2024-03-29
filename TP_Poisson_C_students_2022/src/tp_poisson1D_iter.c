/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <time.h>

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, lab, kv;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *SOL, *EX_SOL, *X;
  double *AB;
  double *MB;
  
  double temp, relres;

  double opt_alpha;

  /* Size of the problem */
  NRHS=1;
  nbpoints=12;
  la=nbpoints-2;

  /* Dirichlet Boundary conditions */
  T0=5.0;
  T1=20.0;

  printf("\n--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  SOL=(double *) calloc(la, sizeof(double)); 
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "./iter dat/RHS.dat");
  write_vec(EX_SOL, &la, "./iter dat/EX_SOL.dat");
  write_vec(X, &la, "./iter dat/X_grid.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  
  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  /* uncomment the following to check matrix A */
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "./iter dat/AB.dat");
  
  /********************************************/
  /* Solution (Richardson with optimal alpha) */
  printf("***** Solution (Richardson with optimal alpha) *****\n");

  /* Computation of optimum alpha */
  opt_alpha = richardson_alpha_opt(&la);
  printf("Optimal alpha for simple Richardson iteration is : %lf\n",opt_alpha); 

  /* Solve */
  double tol=1e-3;
  int maxit=1000;
  double *resvec;
  int nbite=0;

  resvec=(double *) calloc(maxit, sizeof(double));

  /* Solve with Richardson alpha */
  clock_t debut_alpha = clock();
  richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  clock_t fin_alpha = clock();

  printf("Temps d'execution RICHARDSON ALPHA  = %f seconds\n", (double)(fin_alpha - debut_alpha) / CLOCKS_PER_SEC);

  /* Write solution */
  write_vec(SOL, &la, "./iter dat/SOL_ALPHA.dat");


  printf("\n***** l\'erreur RICHARDSON ALPHA par rapport à la solution analytique *****\n");
  /* Relative forward error */
  temp = cblas_ddot(la, SOL, 1, SOL,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, SOL, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;
  
  printf("\nl\'erreur RICHARDSON ALPHA = %f\n",relres);

  
  //liberer pour les utiliser a nouveau
  free(SOL);
  SOL=( double *) calloc(la, sizeof(double));
  nbite = 0; 
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);




  /* Richardson General Tridiag */
  printf("\n***** Richardson General Tridiag *****\n");


  /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
  kv = 1;
  ku = 1;
  kl = 1;
  MB = (double *) malloc(sizeof(double)*(la)*la);

  /* AVEC JACOBI */
  printf("***** Avec JACOBI *****\n");

  extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "./iter dat/MB_JACOBI.dat");

  /* Solve with General Richardson */
  clock_t debut_jacobi = clock();
  richardson_MB_Jacobi(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  clock_t fin_jacobi = clock();

  printf("Temps d'execution RICHARDSON JACOBI  = %f seconds\n", (double)(fin_jacobi - debut_jacobi) / CLOCKS_PER_SEC);
  /* Write solution */
  write_vec(SOL, &la, "./iter dat/SOL_JACOBI.dat");

  printf("\n***** l\'erreur RICHARDSON JACOBI par rapport à la solution analytique *****\n");
  /* Relative forward error */
  temp = cblas_ddot(la, SOL, 1, SOL,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, SOL, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;
  
  printf("\nl\'erreur RICHARDSON JACOBI = %f\n",relres);

  free(SOL);
  SOL=( double *) calloc(la, sizeof(double));
  nbite = 0; 
  free(MB);
  MB = (double *) malloc(sizeof(double)*(la)*la);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);



  /* AVEC GAUSS SEIDEL */
  printf("\n***** Avec GAUSS SEIDEL *****\n");

  extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
  write_GB_operator_colMajor_poisson1D(MB, &la, &la, "./iter dat/MB_GAUSS_SEIDEL.dat");

  /* Solve with General Richardson */
  clock_t debut_gauss = clock();
  richardson_MB_Gauss(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
  clock_t fin_gauss = clock();

  printf("Temps d'execution RICHARDSON GAUSS SEIDEL = %f seconds\n", (double)(fin_gauss - debut_gauss) / CLOCKS_PER_SEC);

  
  /* Write solution */
  write_vec(SOL, &la, "./iter dat/SOL_GAUSS.dat");

  printf("\n***** l\'erreur RICHARDSON GAUSS SEIDEL par rapport à la solution analytique *****\n");
  /* Relative forward error */
  temp = cblas_ddot(la, SOL, 1, SOL,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, SOL, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;
  
  printf("\nl\'erreur RICHARDSON GAUSS SEIDEL = %f\n",relres);

  
  /* Write solution */
  write_vec(SOL, &la, "./iter dat/SOL.dat");

  /* Write convergence history */
  write_vec(resvec, &nbite, "./iter dat/RESVEC.dat");

  free(resvec);
  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  free(MB);
  printf("\n\n--------- End -----------\n");
}
